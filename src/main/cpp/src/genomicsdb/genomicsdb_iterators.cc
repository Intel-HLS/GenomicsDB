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

#define VERIFY_OR_THROW(X) if(!(X)) throw GenomicsDBIteratorException(#X);

SingleCellTileDBIterator::SingleCellTileDBIterator(TileDB_CTX* tiledb_ctx, const VariantArraySchema& variant_array_schema,
        const std::string& array_path, const int64_t* range, const std::vector<int>& attribute_ids, const size_t buffer_size)
: m_variant_array_schema(&variant_array_schema), m_cell(variant_array_schema, attribute_ids),
#ifdef DO_PROFILING
  , m_tiledb_timer()
  , m_tiledb_to_buffer_cell_timer()
#endif
{
  std::vector<const char*> attribute_names(attribute_ids.size()+1u);  //+1 for the COORDS
  m_query_attribute_idx_vec.resize(attribute_ids.size()+1u);//+1 for the COORDS
  m_query_attribute_idx_to_tiledb_buffer_idx.resize(attribute_ids.size()+1u);//+1 for the COORDS
  for(auto i=0ull;i<attribute_ids.size();++i)
  {
    auto schema_idx = attribute_ids[i];
    attribute_names[i] = variant_array_schema.attribute_name(schema_idx).c_str();
    m_query_attribute_idx_vec[i] = i; //initially all attributes are queried
    m_query_attribute_idx_to_tiledb_buffer_idx[i] = m_buffer_pointers.size();
    auto is_variable_length_field = variant_array_schema.is_variable_length_field(schema_idx);
    //GenomicsDBColumnarField
    m_fields.emplace_back(variant_array_schema.type(schema_idx),
        is_variable_length_field ? BCF_VL_VAR : BCF_VL_FIXED,
        variant_array_schema.val_num(schema_idx), buffer_size);
    //Get a free buffer and move it to the live list
    m_fields[i].move_buffer_to_live_list(m_fields[i].get_free_buffer());
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
  m_query_attribute_idx_vec[coords_idx] = coords_idx; //initially all attributes are queried
  m_query_attribute_idx_to_tiledb_buffer_idx[coords_idx] = m_buffer_pointers.size();
  //GenomicsDBColumnarField
  m_fields.emplace_back(variant_array_schema.dim_type(), BCF_VL_FIXED, variant_array_schema.dim_length(), buffer_size);
  m_fields[coords_idx].move_buffer_to_live_list(m_fields[coords_idx].get_free_buffer());
  //Buffer pointers and size for COORDS
  m_buffer_pointers.push_back(0);
  m_buffer_sizes.push_back(0);
  /* Initialize the array in READ mode. */
  auto status = tiledb_array_init(
      tiledb_ctx, 
      &m_tiledb_array,
      array_path.c_str(),
      TILEDB_ARRAY_READ,
      reinterpret_cast<const void*>(range), // range,
      &(attribute_names[0]),           
      attribute_names.size());
  VERIFY_OR_THROW(status == TILEDB_OK && "Error while initializing TileDB array object");
#ifdef DO_PROFILING
  m_tiledb_timer.stop();
#endif
}

void SingleCellTileDBIterator::read_from_TileDB()
{
  //Zero out all buffer sizes
  memset(&(m_buffer_sizes[0]), 0, m_buffer_sizes.size()*sizeof(size_t));
  //Only set non-0 buffer sizes for attributes for which we wish to get more data
  for(auto i=0ull;i<m_query_attribute_idx_vec.size();++i)
  {
    auto query_idx = m_query_attribute_idx_vec[i];
    assert(static_cast<size_t>(query_idx) < m_fields.size());
    auto& genomicsdb_columnar_field = m_fields[query_idx];
    auto* genomicsdb_buffer_ptr = genomicsdb_columnar_field.get_live_buffer_list_head_ptr();
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
  //Set number of live entries in each buffer
  for(auto i=0ull;i<m_query_attribute_idx_vec.size();++i)
  {
    auto query_idx = m_query_attribute_idx_vec[i];
    auto& genomicsdb_columnar_field = m_fields[query_idx];
    auto* genomicsdb_buffer_ptr = genomicsdb_columnar_field.get_live_buffer_list_head_ptr();
    auto buffer_idx = m_query_attribute_idx_to_tiledb_buffer_idx[query_idx];
    //For variable length field, first the offsets buffer
    auto filled_buffer_size = m_buffer_sizes[buffer_idx];
    if(genomicsdb_columnar_field.is_variable_length_field())
    {
      auto num_live_entries = filled_buffer_size/sizeof(size_t);
      genomicsdb_buffer_ptr->set_num_live_entries(num_live_entries);
      //Append size of the variable length field
      genomicsdb_buffer_ptr->set_offset(num_live_entries, m_buffer_sizes[buffer_idx+1u]);
    }
    else
      genomicsdb_buffer_ptr->set_num_live_entries(filled_buffer_size/
            (genomicsdb_columnar_field.get_fixed_length_field_size_in_bytes()));
  }
}
