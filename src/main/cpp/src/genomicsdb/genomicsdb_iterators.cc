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

SingleCellTileDBIterator::SingleCellTileDBIterator(TileDB_CTX* tiledb_ctx,
    const VidMapper* vid_mapper, const VariantArraySchema& variant_array_schema,
    const std::string& array_path, const VariantQueryConfig& query_config, const size_t buffer_size)
  : SingleCellTileDBIterator(tiledb_ctx,
      0,
      vid_mapper, variant_array_schema,
      array_path, query_config, buffer_size)
{}

SingleCellTileDBIterator::SingleCellTileDBIterator(TileDB_CTX* tiledb_ctx,
    const TileDB_Array* tiledb_array,
    const VidMapper* vid_mapper, const VariantArraySchema& variant_array_schema,
    const std::string& array_path, const VariantQueryConfig& query_config, const size_t buffer_size)
: m_variant_array_schema(&variant_array_schema), m_query_config(&query_config),
  m_cell(new GenomicsDBColumnarCell(this)),
  m_tiledb_array(tiledb_array), m_owned_tiledb_array(0),
  m_query_column_interval_idx(0u),
  m_first_read_from_TileDB(true),
  m_done_reading_from_TileDB(false),
  m_in_find_intersecting_intervals_mode(false),
  m_in_simple_traversal_mode(false),
  m_live_cell_markers(query_config.get_num_rows_in_array(), query_config.get_num_queried_attributes()+1u), //+1 for coords
  m_num_markers_initialized(0),
  m_PQ_live_cell_markers(m_live_cell_markers),
  m_smallest_row_idx_in_array(query_config.get_smallest_row_idx_in_array())
#ifdef DO_PROFILING
  , m_tiledb_timer()
  , m_tiledb_to_buffer_cell_timer()
#endif
{
  const auto& attribute_ids = query_config.get_query_attributes_schema_idxs();
  std::vector<const char*> attribute_names(attribute_ids.size()+1u);  //+1 for the COORDS
  m_query_attribute_idx_vec.resize(attribute_ids.size()+1u);//+1 for the COORDS
  m_query_attribute_idx_num_cells_to_increment_vec.resize(attribute_ids.size()+1u, 0u);
  m_query_attribute_idx_to_tiledb_buffer_idx.resize(attribute_ids.size()+1u);//+1 for the COORDS
  for(auto i=0ull;i<attribute_ids.size();++i)
  {
    auto schema_idx = attribute_ids[i];
    attribute_names[i] = variant_array_schema.attribute_name(schema_idx).c_str();
    m_query_attribute_idx_to_tiledb_buffer_idx[i] = m_buffer_pointers.size();
    auto is_variable_length_field = variant_array_schema.is_variable_length_field(schema_idx);
    //TileDB does not have a way to distinguish between char, string and int8_t fields
    //Hence, a series of possibly messy checks here
    assert(vid_mapper);
    const auto* vid_field_info = vid_mapper->get_field_info(attribute_names[i]);
    auto curr_type_index = (vid_field_info && vid_field_info->m_bcf_ht_type == BCF_HT_FLAG)
      ? std::type_index(typeid(bool))
      : variant_array_schema.type(schema_idx);
    //GenomicsDBColumnarField
    m_fields.emplace_back(curr_type_index,
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
  //Set row idx for live cell markers
  for(auto i=0ull;i<query_config.get_num_rows_in_array();++i)
    m_live_cell_markers.set_row_idx(i, m_smallest_row_idx_in_array+i);
  //END query idx
  assert(m_query_config->is_defined_query_idx_for_known_field_enum(GVCF_END_IDX));
  m_END_query_idx = m_query_config->get_query_idx_for_known_field_enum(GVCF_END_IDX);
  assert(m_END_query_idx < m_fields.size() && m_END_query_idx == 0u); //END must always be the first field
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
  if(m_owned_tiledb_array)
    tiledb_array_finalize(m_owned_tiledb_array);
  m_owned_tiledb_array = 0;
  m_tiledb_array = 0;
#ifdef DO_PROFILING
  m_tiledb_timer.print("TileDB iterator", std::cerr);
  m_tiledb_to_buffer_cell_timer.print("TileDB to buffer cell", std::cerr);
#endif
}

void SingleCellTileDBIterator::read_from_TileDB(TileDB_CTX* tiledb_ctx, const char* array_path,
    std::vector<const char*>* attribute_names)
{
  //This function gets called recursively
  auto this_stack_frame_in_find_intersecting_intervals_mode = false;
  do
  {
    //First time or completed query for m_query_column_interval_idx or
    //this stack was in find intersecting mode in the previous iteration
    if(m_first_read_from_TileDB || m_done_reading_from_TileDB || this_stack_frame_in_find_intersecting_intervals_mode)
    {
      //shouldn't enter this if condition if called recursively through operator++ in m_in_find_intersecting_intervals_mode
      assert(!m_in_find_intersecting_intervals_mode);
      //this stack was in find intersecting mode in the previous iteration
      //Move into simple traversal mode
      if(this_stack_frame_in_find_intersecting_intervals_mode)
      {
        m_done_reading_from_TileDB = false;
        this_stack_frame_in_find_intersecting_intervals_mode = false;
        m_in_simple_traversal_mode = true;
      }
      if(m_done_reading_from_TileDB) //move to the next query interval
      {
        m_done_reading_from_TileDB = false;
        m_in_simple_traversal_mode = false;
        ++m_query_column_interval_idx;
        assert(m_query_column_interval_idx < m_query_config->get_num_column_intervals());
        //Reset markers
        assert(m_PQ_live_cell_markers.empty());
        m_live_cell_markers.reset();
        m_num_markers_initialized = 0u;
        //Move all buffers to free list - moving to the next query interval
        for(auto& field : m_fields)
          field.move_all_buffers_from_live_list_to_free_list();
      }
      //Query all attributes
      m_query_attribute_idx_vec.resize(m_fields.size());
      for(auto i=0u;i<m_fields.size();++i)
        m_query_attribute_idx_vec[i] = i;
      //Range of the query
      int64_t query_range[4] = { m_smallest_row_idx_in_array,
        static_cast<int64_t>(m_query_config->get_num_rows_in_array()+m_smallest_row_idx_in_array-1), //inclusive bounds
        0, INT64_MAX-1 };
      if(m_query_config->get_num_column_intervals() > 0u)
      {
        query_range[2] = m_query_config->get_column_begin(m_query_column_interval_idx);
        //Done with find intersecting mode, now just simple traversal
        if(m_in_simple_traversal_mode)
          query_range[3] = m_query_config->get_column_end(m_query_column_interval_idx);
        else
        {
          this_stack_frame_in_find_intersecting_intervals_mode = true;
          m_in_find_intersecting_intervals_mode = true;
        }
      }
      else
        m_in_simple_traversal_mode = true; //no query intervals, full array scan, simple traversal
      auto status = -1;
      if(m_tiledb_array == 0)
      {
        assert(array_path && attribute_names && tiledb_ctx);
        /* Initialize the array in READ mode. */
        status = tiledb_array_init(
            tiledb_ctx,
            &m_owned_tiledb_array,
            array_path,
            TILEDB_ARRAY_READ,
            reinterpret_cast<const void*>(query_range),
            &((*attribute_names)[0]),
            attribute_names->size());
        VERIFY_OR_THROW(status == TILEDB_OK && "Error while initializing TileDB array object");
        m_tiledb_array = m_owned_tiledb_array;
      }
      else //next column interval - reset subarray
      {
        //TileDB_Array object is provided by caller and this is the first read, reset attributes to query
        if(m_first_read_from_TileDB)
        {
          status = tiledb_array_reset_attributes(m_tiledb_array,
              &((*attribute_names)[0]),
              attribute_names->size());
          VERIFY_OR_THROW(status == TILEDB_OK && "Error while initializing attributes for the TileDB array object");
        }
        status = tiledb_array_reset_subarray(m_tiledb_array, reinterpret_cast<const void*>(query_range));
        VERIFY_OR_THROW(status == TILEDB_OK && "Error in tiledb_array_reset_subarray()");
      }
      m_first_read_from_TileDB = false;
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
      genomicsdb_columnar_field.move_buffer_to_live_list(genomicsdb_columnar_field.get_or_allocate_free_buffer());
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
        genomicsdb_columnar_field.move_buffer_to_free_list(genomicsdb_buffer_ptr);
      }
    }
#ifdef DEBUG
    //The way the code is written, the final call to read_from_TileDB() will attempt to
    //fetch data for all queried fields. TileDB should return 0 for all fields. Hence, either all
    //queried fields return with 0 size or have some data to process
    assert(num_done_fields == 0u || (num_done_fields == m_query_attribute_idx_vec.size() && m_done_reading_from_TileDB));
#endif
    if(this_stack_frame_in_find_intersecting_intervals_mode)
    {
      while(!m_done_reading_from_TileDB)
      {
        handle_current_cell_in_find_intersecting_intervals_mode();
        if(m_in_find_intersecting_intervals_mode)
          this->operator++();
        else
          break;
      }
      m_in_find_intersecting_intervals_mode = false;
    }
    //If done reading current column interval, move to next query column interval, if any remaining OR
    //this stack was performing intersecting intervals search
  } while((m_done_reading_from_TileDB && m_PQ_live_cell_markers.empty()
        && (m_query_column_interval_idx+1u) < m_query_config->get_num_column_intervals())
      || this_stack_frame_in_find_intersecting_intervals_mode);
}

void SingleCellTileDBIterator::handle_current_cell_in_find_intersecting_intervals_mode()
{
  assert(m_in_find_intersecting_intervals_mode && !m_done_reading_from_TileDB);
  auto& coords_columnar_field = m_fields[m_fields.size()-1u];
  auto* coords = reinterpret_cast<int64_t*>(
      coords_columnar_field.get_pointer_to_data_in_buffer_at_index(
        coords_columnar_field.get_live_buffer_list_tail_ptr(),
        coords_columnar_field.get_curr_index_in_live_list_tail()
        )
      );
  auto& END_columnar_field = m_fields[m_END_query_idx];
  auto& END_field_value = *(reinterpret_cast<int64_t*>(
        END_columnar_field.get_pointer_to_data_in_buffer_at_index(
          END_columnar_field.get_live_buffer_list_tail_ptr(),
          END_columnar_field.get_curr_index_in_live_list_tail()
          )
        ));
  auto row_idx = coords[0];
  assert(row_idx >= m_smallest_row_idx_in_array);
  auto marker_idx = row_idx - m_smallest_row_idx_in_array;
  assert(m_live_cell_markers.get_row_idx(marker_idx) == row_idx);
  //Not yet initialized
  if(!m_live_cell_markers.is_initialized(marker_idx))
  {
    m_live_cell_markers.set_initialized(marker_idx, true); //one way or the other this is initialized
    auto& coords_column = coords[1];
    assert(m_query_column_interval_idx < m_query_config->get_num_column_intervals());
    auto query_interval_begin = m_query_config->get_column_begin(m_query_column_interval_idx);
    //Only deal with END cells that begin before the queried column interval begin and intersect it
    if(coords_column > END_field_value && static_cast<uint64_t>(END_field_value) < query_interval_begin)
    {
      std::swap(coords_column, END_field_value);
      m_live_cell_markers.set_valid(marker_idx, true);
      //Add to PQ
      m_live_cell_markers.set_column_interval(marker_idx, coords_column, END_field_value);
      m_PQ_live_cell_markers.push(marker_idx);
      //Track offsets for each field and keep the buffer alive
      for(auto i=0u;i<m_fields.size();++i)
      {
        auto& genomicsdb_columnar_field = m_fields[i];
        auto* genomicsdb_buffer_ptr = genomicsdb_columnar_field.get_live_buffer_list_tail_ptr();
        assert(genomicsdb_buffer_ptr);
        m_live_cell_markers.set_field_marker(marker_idx, i, genomicsdb_buffer_ptr,
            genomicsdb_columnar_field.get_curr_index_in_live_list_tail());
        genomicsdb_buffer_ptr->increment_num_live_entries();
      }
    }
    else
      m_live_cell_markers.set_valid(marker_idx, false);
    ++m_num_markers_initialized;
    //All initialized
    if(m_num_markers_initialized == m_query_config->get_num_rows_in_array())
      m_in_find_intersecting_intervals_mode = false;
  }
}

const SingleCellTileDBIterator& SingleCellTileDBIterator::operator++()
{
  auto curr_query_column_interval_idx = m_query_column_interval_idx;
  //In simple traversal mode, but intervals intersecting the column begin still exist in PQ
  if(m_in_simple_traversal_mode && !m_PQ_live_cell_markers.empty())
  {
    auto marker_idx = m_PQ_live_cell_markers.top();
    m_PQ_live_cell_markers.pop();
    //Decrease #live entries and move to free list if needed
    for(auto i=0u;i<m_fields.size();++i)
    {
      auto* buffer_ptr = m_live_cell_markers.get_buffer_pointer(marker_idx, i);
      buffer_ptr->decrement_num_live_entries(1u);
      if(buffer_ptr->get_num_live_entries() == 0u)
        m_fields[i].move_buffer_to_free_list(buffer_ptr);
    }
    if(!m_PQ_live_cell_markers.empty())
      return *this;
    //Done processing PQ cells
    //If this query interval has no more data, read next column interval
    if(m_done_reading_from_TileDB
        && (curr_query_column_interval_idx+1u) < m_query_config->get_num_column_intervals())
    {
      read_from_TileDB();
      return *this;
    }
  }
  auto hitting_useless_cells = true;
  //Must increment by 1 at least
  auto num_cells_incremented = 1ull;
  const auto& coords_columnar_field = m_fields[m_fields.size()-1u];
  assert(m_END_query_idx < m_fields.size());
  const auto& END_columnar_field = m_fields[m_END_query_idx];
  while(hitting_useless_cells && !m_done_reading_from_TileDB)
  {
    //Increment coords and END fields only first - increment other fields only when a useful cell is found
    m_query_attribute_idx_vec.resize(2u); //no heap operations here
    m_query_attribute_idx_num_cells_to_increment_vec.resize(2u);
    m_query_attribute_idx_vec[0u] = m_fields.size()-1u; //coords
    m_query_attribute_idx_vec[1u] = m_END_query_idx;
    m_query_attribute_idx_num_cells_to_increment_vec[0u] = 1u;  //single cell increment
    m_query_attribute_idx_num_cells_to_increment_vec[1u] = 1u;
    increment_iterator_within_live_buffer_list_tail_ptr_for_fields();
    if(m_query_attribute_idx_vec.size() > 0u) //some fields have exhausted buffers, need to fetch from TileDB
    {
      read_from_TileDB();
      //The read_from_TileDB() might have determined that no more cells exist for the current query
      //interval. It would move on to the next interval. So, this stack frame must assume the query is done
      //and do nothing
      if(m_done_reading_from_TileDB
          || (curr_query_column_interval_idx != m_query_column_interval_idx))
        return *this;
    }
    const auto* coords = reinterpret_cast<const int64_t*>(
        coords_columnar_field.get_pointer_to_data_in_buffer_at_index(
          coords_columnar_field.get_live_buffer_list_tail_ptr(),
          coords_columnar_field.get_curr_index_in_live_list_tail()
          )
        );
    //keep incrementing iterator if hitting duplicates at end in simple traversal mode
    if(m_in_simple_traversal_mode)
    {
      assert(END_columnar_field.get_live_buffer_list_tail_ptr()->get_num_unprocessed_entries() > 0u);
      auto END_field_value = *(reinterpret_cast<const int64_t*>(
            END_columnar_field.get_pointer_to_data_in_buffer_at_index(
              END_columnar_field.get_live_buffer_list_tail_ptr(),
              END_columnar_field.get_curr_index_in_live_list_tail()
              )
            ));
      hitting_useless_cells = (END_field_value < coords[1]);
    }
    else
    {
      assert(m_in_find_intersecting_intervals_mode);
      auto row_idx = coords[0];
      assert(row_idx >= m_smallest_row_idx_in_array);
      auto marker_idx = row_idx - m_smallest_row_idx_in_array;
      assert(m_live_cell_markers.get_row_idx(marker_idx) == row_idx);
      //if the row is already initialized in the find intersecting intervals mode, can skip
      //this cell
      hitting_useless_cells = m_live_cell_markers.is_initialized(marker_idx);
    }
    if(hitting_useless_cells)
      ++num_cells_incremented;
  }
  assert(m_in_find_intersecting_intervals_mode || m_PQ_live_cell_markers.empty());
  //Increment iterator for other fields
  //All fields other than coords and END
  m_query_attribute_idx_vec.resize(m_fields.size()-2u);
  //END must be the first field
  assert(m_END_query_idx == 0u);
  for(auto i=0u;i<m_query_attribute_idx_vec.size();++i)
    m_query_attribute_idx_vec[i] = i+1u;  //END is the first field - ignore
  //For all fields, #cells to skip == num_cells_incremented initially
  m_query_attribute_idx_num_cells_to_increment_vec.resize(m_fields.size()-2u);
  m_query_attribute_idx_num_cells_to_increment_vec.assign(
      m_query_attribute_idx_num_cells_to_increment_vec.size(), num_cells_incremented);
  increment_iterator_within_live_buffer_list_tail_ptr_for_fields();
  //some fields have exhausted buffers, need to fetch from TileDB
  while(!m_query_attribute_idx_vec.empty())
  {
    read_from_TileDB();
    //TODO: either all fields are done reading or none are - is that right?
    //also, same column interval being queried
    assert(!m_done_reading_from_TileDB
        && m_query_column_interval_idx == curr_query_column_interval_idx);
    increment_iterator_within_live_buffer_list_tail_ptr_for_fields();
  }
#ifdef DEBUG
  for(auto i=0u;i<m_fields.size();++i)
    assert(m_fields[i].get_live_buffer_list_tail_ptr()
        && m_fields[i].get_live_buffer_list_tail_ptr()->get_num_unprocessed_entries());
#endif
  return *this;
}

void SingleCellTileDBIterator::print(const int field_query_idx, std::ostream& fptr) const
{
  assert(static_cast<const size_t>(field_query_idx) < m_fields.size());
  auto& genomicsdb_columnar_field = m_fields[field_query_idx];
  size_t index = 0ul;
  auto buffer_ptr = get_buffer_pointer_and_index(field_query_idx, index);
  genomicsdb_columnar_field.print_data_in_buffer_at_index(fptr,
      buffer_ptr, index);
}

void SingleCellTileDBIterator::increment_iterator_within_live_buffer_list_tail_ptr_for_fields()
{
  auto next_iteration_num_query_attributes = 0u;
  for(auto i=0u;i<m_query_attribute_idx_vec.size();++i)
  {
    assert(next_iteration_num_query_attributes <= i);
    assert(static_cast<size_t>(m_query_attribute_idx_vec[i]) < m_fields.size());
    auto& genomicsdb_columnar_field = m_fields[m_query_attribute_idx_vec[i]];
    auto* genomicsdb_buffer_ptr = genomicsdb_columnar_field.get_live_buffer_list_tail_ptr();
    //TODO: either all fields are done reading or none are - is that right?
    if(genomicsdb_buffer_ptr)
    {
      auto unprocessed_entries_delta = std::min(m_query_attribute_idx_num_cells_to_increment_vec[i],
          genomicsdb_buffer_ptr->get_num_unprocessed_entries());
      genomicsdb_columnar_field.advance_curr_index_in_live_list_tail(unprocessed_entries_delta);
      genomicsdb_buffer_ptr->decrement_num_live_entries(unprocessed_entries_delta);
      genomicsdb_buffer_ptr->decrement_num_unprocessed_entries(unprocessed_entries_delta);
      if(genomicsdb_buffer_ptr->get_num_live_entries() == 0ull) //no more live entries, move buffer to free list 
        genomicsdb_columnar_field.move_buffer_to_free_list(genomicsdb_buffer_ptr);
      //buffer completely processed, add attribute to next round of fetch from TileDB
      if(genomicsdb_buffer_ptr->get_num_unprocessed_entries() == 0ull)
      {
        m_query_attribute_idx_vec[next_iteration_num_query_attributes]
          = m_query_attribute_idx_vec[i];
        m_query_attribute_idx_num_cells_to_increment_vec[next_iteration_num_query_attributes]
          = m_query_attribute_idx_num_cells_to_increment_vec[i] - unprocessed_entries_delta;
        ++next_iteration_num_query_attributes;
      }
    }
  }
  m_query_attribute_idx_vec.resize(next_iteration_num_query_attributes); //no heap operations occur here
}

void SingleCellTileDBIterator::print_ALT(const int field_query_idx, std::ostream& fptr) const
{
  assert(static_cast<const size_t>(field_query_idx) < m_fields.size());
  auto& genomicsdb_columnar_field = m_fields[field_query_idx];
  size_t index = 0ul;
  auto buffer_ptr = get_buffer_pointer_and_index(field_query_idx, index);
  genomicsdb_columnar_field.print_ALT_data_in_buffer_at_index(fptr,
      buffer_ptr, index);
}

void SingleCellTileDBIterator::print_csv(const int field_query_idx, std::ostream& fptr) const
{
  assert(static_cast<const size_t>(field_query_idx) < m_fields.size());
  auto& genomicsdb_columnar_field = m_fields[field_query_idx];
  size_t index = 0ul;
  auto buffer_ptr = get_buffer_pointer_and_index(field_query_idx, index);
  genomicsdb_columnar_field.print_data_in_buffer_at_index_as_csv(fptr,
      buffer_ptr, index);
}
