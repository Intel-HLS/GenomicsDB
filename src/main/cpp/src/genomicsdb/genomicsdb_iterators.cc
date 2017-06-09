/**
 * The MIT License (MIT)
 * Copyright (c) 2016 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of 
 * this software and associated documentation files (the "Software"), to deal in 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
 * the Software, and to permit persons to whom the Software is furnished to do so, 
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all 
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "genomicsdb_iterators.h"
#include "variant_cell.h"
#include "variant_query_config.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw GenomicsDBIteratorException(#X);

SingleCellTileDBIterator::SingleCellTileDBIterator(TileDB_CTX* tiledb_ctx, const VariantArraySchema& variant_array_schema,
        const std::string& array_path, const VariantQueryConfig& query_config, const size_t buffer_size)
: m_variant_array_schema(&variant_array_schema), m_query_config(&query_config),
  m_done_reading_from_TileDB(false),
  m_cell(new GenomicsDBColumnarCell(this)),
  m_query_column_interval_idx(0u),
  m_tiledb_array(0)
#ifdef DO_PROFILING
  , m_tiledb_timer()
  , m_tiledb_to_buffer_cell_timer()
#endif
{
  const auto& attribute_ids = query_config.get_query_attributes_schema_idxs();
  std::vector<const char*> attribute_names(attribute_ids.size()+1u);  //+1 for the COORDS
  m_query_attribute_idx_vec.resize(attribute_ids.size()+1u);//+1 for the COORDS
  m_query_attribute_idx_to_tiledb_buffer_idx.resize(attribute_ids.size()+1u);//+1 for the COORDS
  for(auto i=0ull;i<attribute_ids.size();++i)
  {
    auto schema_idx = attribute_ids[i];
    attribute_names[i] = variant_array_schema.attribute_name(schema_idx).c_str();
    m_query_attribute_idx_to_tiledb_buffer_idx[i] = m_buffer_pointers.size();
    auto is_variable_length_field = variant_array_schema.is_variable_length_field(schema_idx);
    //GenomicsDBColumnarField
    m_fields.emplace_back(variant_array_schema.type(schema_idx),
        is_variable_length_field ? BCF_VL_VAR : BCF_VL_FIXED,
        variant_array_schema.val_num(schema_idx), buffer_size);
    //Buffer pointers and size
    m_buffer_pointers.push_back(0);
    m_buffer_sizes.push_back(0);
    if(is_variable_length_field)
    {
      m_buffer_pointers.push_back(0);
      m_buffer_sizes.push_back(0);
    }
  }
  //Co-ordinates
  auto coords_idx = attribute_ids.size();
  attribute_names[coords_idx] = TILEDB_COORDS;
  m_query_attribute_idx_to_tiledb_buffer_idx[coords_idx] = m_buffer_pointers.size();
  //GenomicsDBColumnarField
  m_fields.emplace_back(variant_array_schema.dim_type(), BCF_VL_FIXED, variant_array_schema.dim_length(), buffer_size);
  //Buffer pointers and size for COORDS
  m_buffer_pointers.push_back(0);
  m_buffer_sizes.push_back(0);
  //first read
  read_from_TileDB(tiledb_ctx, array_path.c_str(), &attribute_names);
#ifdef DO_PROFILING
  m_tiledb_timer.stop();
#endif
}

SingleCellTileDBIterator::~SingleCellTileDBIterator()
{
  if(m_cell)
    delete m_cell;
  m_cell = 0;
  if(m_tiledb_array)
    tiledb_array_finalize(m_tiledb_array);
  m_tiledb_array = 0;
#ifdef DO_PROFILING
  m_tiledb_timer.print("TileDB iterator", std::cerr);
  m_tiledb_to_buffer_cell_timer.print("TileDB to buffer cell", std::cerr);
#endif
}

void SingleCellTileDBIterator::read_from_TileDB(TileDB_CTX* tiledb_ctx, const char* array_path,
    std::vector<const char*>* attribute_names)
{
  do
  {
    //First time or completed query for m_query_column_interval_idx
    if(m_tiledb_array == 0 || m_done_reading_from_TileDB)
    {
      if(m_done_reading_from_TileDB) //move to the next query interval
      {
        m_done_reading_from_TileDB = false;
        ++m_query_column_interval_idx;
        assert(m_query_column_interval_idx < m_query_config->get_num_column_intervals());
      }
      //Query all attributes
      m_query_attribute_idx_vec.resize(m_fields.size());
      for(auto i=0u;i<m_fields.size();++i)
        m_query_attribute_idx_vec[i] = i;
      //Range of the query
      int64_t query_range[4] = { m_query_config->get_smallest_row_idx_in_array(),
        static_cast<int64_t>(m_query_config->get_num_rows_in_array()+m_query_config->get_smallest_row_idx_in_array()-1),
        0, INT64_MAX-1 };
      if(m_query_config->get_num_column_intervals() > 0u)
      {
        query_range[2] = m_query_config->get_column_begin(m_query_column_interval_idx);
        query_range[3] = m_query_config->get_column_end(m_query_column_interval_idx);
      }
      auto status = -1;
      if(m_tiledb_array == 0)
      {
        assert(array_path && attribute_names && tiledb_ctx);
        /* Initialize the array in READ mode. */
        status = tiledb_array_init(
            tiledb_ctx,
            &m_tiledb_array,
            array_path,
            TILEDB_ARRAY_READ,
            reinterpret_cast<const void*>(query_range),
            &((*attribute_names)[0]),
            attribute_names->size());
        VERIFY_OR_THROW(status == TILEDB_OK && "Error while initializing TileDB array object");
      }
      else //next column interval - reset subarray
      {
        status = tiledb_array_reset_subarray(m_tiledb_array, reinterpret_cast<const void*>(query_range));
        VERIFY_OR_THROW(status == TILEDB_OK && "Error in tiledb_array_reset_subarray()");
      }
    }
    //Zero out all buffer sizes
    memset(&(m_buffer_sizes[0]), 0, m_buffer_sizes.size()*sizeof(size_t));
    //Only set non-0 buffer sizes for attributes for which we wish to get more data
    for(auto i=0ull;i<m_query_attribute_idx_vec.size();++i)
    {
      auto query_idx = m_query_attribute_idx_vec[i];
      assert(static_cast<size_t>(query_idx) < m_fields.size());
      auto& genomicsdb_columnar_field = m_fields[query_idx];
      //Get a free buffer and move it to the live list - buffer gets added to tail of live list
      genomicsdb_columnar_field.move_buffer_to_live_list(genomicsdb_columnar_field.get_free_buffer());
      //Always read new data into tail
      auto* genomicsdb_buffer_ptr = genomicsdb_columnar_field.get_live_buffer_list_tail_ptr();
      assert(genomicsdb_buffer_ptr);
      assert(static_cast<size_t>(query_idx) < m_query_attribute_idx_to_tiledb_buffer_idx.size());
      auto buffer_idx = m_query_attribute_idx_to_tiledb_buffer_idx[query_idx];
      assert(buffer_idx < m_buffer_pointers.size());
      //For variable length field, first the offsets buffer
      if(genomicsdb_columnar_field.is_variable_length_field())
      {
        m_buffer_pointers[buffer_idx] = reinterpret_cast<void*>(genomicsdb_buffer_ptr->get_offsets_pointer());
        m_buffer_sizes[buffer_idx] = genomicsdb_buffer_ptr->get_offsets_size_in_bytes();
        ++buffer_idx;
        assert(buffer_idx < m_buffer_pointers.size());
      }
      m_buffer_pointers[buffer_idx] = reinterpret_cast<void*>(genomicsdb_buffer_ptr->get_buffer_pointer());
      m_buffer_sizes[buffer_idx] = genomicsdb_buffer_ptr->get_buffer_size_in_bytes();
    }
    auto status = tiledb_array_read(m_tiledb_array, &(m_buffer_pointers[0]), &(m_buffer_sizes[0]));
    VERIFY_OR_THROW(status == TILEDB_OK);
#ifdef DEBUG
    auto num_done_fields = 0u;
#endif
    //Set number of live entries in each buffer
    for(auto i=0ull;i<m_query_attribute_idx_vec.size();++i)
    {
      auto query_idx = m_query_attribute_idx_vec[i];
      auto& genomicsdb_columnar_field = m_fields[query_idx];
      auto* genomicsdb_buffer_ptr = genomicsdb_columnar_field.get_live_buffer_list_tail_ptr();
      auto buffer_idx = m_query_attribute_idx_to_tiledb_buffer_idx[query_idx];
      //For variable length field, first the offsets buffer
      auto filled_buffer_size = m_buffer_sizes[buffer_idx];
      auto num_live_entries = 0u;
      if(genomicsdb_columnar_field.is_variable_length_field())
      {
        num_live_entries = filled_buffer_size/sizeof(size_t);
        genomicsdb_buffer_ptr->set_num_live_entries(num_live_entries);
        //Append size of the variable length field
        genomicsdb_buffer_ptr->set_offset(num_live_entries, m_buffer_sizes[buffer_idx+1u]);
      }
      else
      {
        num_live_entries = filled_buffer_size/
          (genomicsdb_columnar_field.get_fixed_length_field_size_in_bytes());
        genomicsdb_buffer_ptr->set_num_live_entries(num_live_entries);
      }
      //Marks data as valid/invalid
      genomicsdb_columnar_field.set_valid_vector_in_live_buffer_list_tail_ptr();
      if(num_live_entries == 0u)
      {
#ifdef DEBUG
        ++num_done_fields;
#endif
        m_done_reading_from_TileDB = true;
      }
    }
#ifdef DEBUG
    //The way the code is written, the final call to read_from_TileDB() will attempt to
    //fetch data for all fields. TileDB should return 0 for all fields. Hence, either all
    //fields are queried and return with 0 size or none do
    assert(num_done_fields == 0u || (num_done_fields == m_fields.size() && m_done_reading_from_TileDB));
#endif
    //If done reading current column interval, move to next query column interval, if any remaining
  } while(m_done_reading_from_TileDB && (m_query_column_interval_idx+1u) < m_query_config->get_num_column_intervals());
}

const SingleCellTileDBIterator& SingleCellTileDBIterator::operator++()
{
  assert(m_query_config->is_defined_query_idx_for_known_field_enum(GVCF_END_IDX));
  auto END_query_idx = m_query_config->get_query_idx_for_known_field_enum(GVCF_END_IDX);
  assert(END_query_idx < m_fields.size());
  auto hitting_duplicate_cells_at_end = true;
  while(hitting_duplicate_cells_at_end && !m_done_reading_from_TileDB)
  {
    m_query_attribute_idx_vec.resize(m_fields.size()); //no heap operations here
    auto next_iteration_num_query_attributes = 0u;
    for(auto i=0u;i<m_fields.size();++i)
    {
      auto& genomicsdb_columnar_field = m_fields[i];
      auto* genomicsdb_buffer_ptr = genomicsdb_columnar_field.get_live_buffer_list_tail_ptr();
      //still has live data
      //TODO: either all fields have live data or none do - is that right?
      if(genomicsdb_buffer_ptr)
      {
        genomicsdb_columnar_field.advance_curr_index_in_live_list_tail();
        genomicsdb_buffer_ptr->decrement_num_live_entries();
        if(genomicsdb_buffer_ptr->get_num_live_entries() == 0ull) //no more live entries, add to queries in next round
        {
          genomicsdb_columnar_field.move_buffer_to_free_list(genomicsdb_buffer_ptr);
          m_query_attribute_idx_vec[next_iteration_num_query_attributes++] = i;
        }
      }
    }
    m_query_attribute_idx_vec.resize(next_iteration_num_query_attributes); //no heap operations occur here
    if(next_iteration_num_query_attributes > 0u) //some fields have exhausted buffers, need to fetch from TileDB
      read_from_TileDB();
    //keep incrementing iterator if hitting duplicates at end
    if(!m_done_reading_from_TileDB)
    {
      const auto& coords_columnar_field = m_fields[m_fields.size()-1u];
      const auto* coords = reinterpret_cast<const int64_t*>(
          coords_columnar_field.get_pointer_to_data_in_buffer_at_index(
            coords_columnar_field.get_live_buffer_list_tail_ptr(),
            coords_columnar_field.get_curr_index_in_live_list_tail()
            )
          );
      const auto& END_columnar_field = m_fields[END_query_idx];
      auto END_position = *(reinterpret_cast<const int64_t*>(
            END_columnar_field.get_pointer_to_data_in_buffer_at_index(
              END_columnar_field.get_live_buffer_list_tail_ptr(),
              END_columnar_field.get_curr_index_in_live_list_tail()
              )
            ));
      if(END_position >= coords[1])
        hitting_duplicate_cells_at_end = false;
    }
  }
  return *this;
}

void SingleCellTileDBIterator::print(const int query_idx, std::ostream& fptr) const
{
  assert(static_cast<const size_t>(query_idx) < m_fields.size());
  auto& genomicsdb_columnar_field = m_fields[query_idx];
  genomicsdb_columnar_field.print_data_in_buffer_at_index(fptr,
      genomicsdb_columnar_field.get_live_buffer_list_tail_ptr(),
      genomicsdb_columnar_field.get_curr_index_in_live_list_tail());
}
