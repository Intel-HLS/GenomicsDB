/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

#include "variant_cell.h"
#include "variant_field_data.h"
#include "variant_query_config.h"
#include "json_config.h"

BufferVariantCell::BufferVariantCell(const VariantArraySchema& array_schema)
{
  clear();
  m_array_schema = &array_schema;
  resize(array_schema.attribute_num());
  for(auto i=0u;i<array_schema.attribute_num();++i)
    update_field_info(i, i);
}

BufferVariantCell::BufferVariantCell(const VariantArraySchema& array_schema, const VariantQueryConfig& query_config)
{
  clear();
  m_array_schema = &array_schema;
  assert(query_config.is_bookkeeping_done());
  resize(query_config.get_num_queried_attributes());
  for(auto i=0u;i<query_config.get_num_queried_attributes();++i)
    update_field_info(i, query_config.get_schema_idx_for_query_idx(i));
}

BufferVariantCell::BufferVariantCell(const VariantArraySchema& array_schema, const std::vector<int>& attribute_ids)
{
  clear();
  m_array_schema = &array_schema;
  resize(attribute_ids.size());
  for(auto i=0u;i<attribute_ids.size();++i)
    update_field_info(i, attribute_ids[i]);
}

void BufferVariantCell::clear()
{
  m_schema_idxs.clear();
  m_field_ptrs.clear();
  m_field_lengths.clear();
  m_row_idx = m_begin_column_idx = -1ll;
}

void BufferVariantCell::resize(const size_t num_fields)
{
  //Resize vectors
  m_schema_idxs.resize(num_fields);
  m_field_ptrs.resize(num_fields);
  m_field_lengths.resize(num_fields);
}
    
void BufferVariantCell::update_field_info(const int query_idx, const int schema_idx)
{
  assert(static_cast<const size_t>(query_idx) < m_field_ptrs.size() && m_field_ptrs.size() == m_schema_idxs.size()
      && m_field_lengths.size() == m_schema_idxs.size());
  m_schema_idxs[query_idx] = schema_idx;
  m_field_lengths[query_idx] = m_array_schema->val_num(schema_idx);
}

void BufferVariantCell::set_cell(const void* ptr)
{
  assert(ptr);
  auto cell_ptr = reinterpret_cast<const uint8_t*>(ptr);
  m_row_idx = *(reinterpret_cast<const int64_t*>(cell_ptr));
  m_begin_column_idx = *(reinterpret_cast<const int64_t*>(cell_ptr+sizeof(int64_t))); 
  //Go past co-ordinates and cell size
  uint64_t offset = 2*sizeof(int64_t)+sizeof(size_t);
#ifdef DEBUG
  auto cell_size = *(reinterpret_cast<const size_t*>(cell_ptr+2*sizeof(int64_t)));
#endif
  for(auto i=0u;i<m_field_ptrs.size();++i)
  {
    auto schema_idx = m_schema_idxs[i];
    auto length = m_field_lengths[i];
    //#define DEBUG_VARIANT_CELL_OFFSETS
#ifdef DEBUG_VARIANT_CELL_OFFSETS
    //For variable length fields, data_offset == length_offset+sizeof(int)
    //For fixed size fields, data_offset == length_offset
    std::cerr << m_array_schema->attribute_name(schema_idx)
      << " length_offset "<< offset;
#endif
    //check if variable length field - read length from buffer
    if(m_array_schema->is_variable_length_field(schema_idx))
    {
      length = *(reinterpret_cast<const int*>(cell_ptr+offset));
      m_field_lengths[i] = length;
      offset += sizeof(int);
    }
#ifdef DEBUG_VARIANT_CELL_OFFSETS
    std::cerr << " data_offset "<<offset<<" #elements "<< length << "\n";
#endif
    m_field_ptrs[i] = cell_ptr + offset;      //field pointer points to region in buffer AFTER the length
    offset += (length*(m_array_schema->element_size(schema_idx)));
  }
#ifdef DEBUG
  assert(offset == cell_size);
#endif
}

void GenomicsDBColumnarCell::print(std::ostream& fptr, const VariantQueryConfig* query_config,
    const std::string& indent_prefix, const VidMapper* vid_mapper) const
{
  assert(query_config && query_config->is_bookkeeping_done());
  auto coords = get_coordinates();
  auto indent_string = indent_prefix+g_json_indent_unit;
  fptr << indent_prefix << "{\n";
  fptr << indent_string << "\"row\": "<< coords[0] << ",\n";
  assert(query_config->is_defined_query_idx_for_known_field_enum(GVCF_END_IDX));
  auto end_query_idx = query_config->get_query_idx_for_known_field_enum(GVCF_END_IDX);
  auto ALT_query_idx = query_config->is_defined_query_idx_for_known_field_enum(GVCF_ALT_IDX)
    ? query_config->get_query_idx_for_known_field_enum(GVCF_ALT_IDX) : UINT64_MAX;
  auto end_position = *(reinterpret_cast<const int64_t*>(get_field_ptr_for_query_idx(end_query_idx)));
  fptr << indent_string << "\"interval\": [ "<< coords[1] << ", "
    << end_position 
    << " ],\n";
  assert(coords[1] <= end_position);
  if(vid_mapper)
  {
    std::string contig_name;
    int64_t contig_position;
    auto status = vid_mapper->get_contig_location(coords[1], contig_name, contig_position);
    if(status)
      fptr << indent_string << "\"genomic_interval\": { \"" << contig_name << "\" : [ "<<contig_position+1
        << ", " << (contig_position+1+(end_position - coords[1])) << " ] },\n";
  }
  fptr << indent_string << "\"fields\": {\n";
  indent_string += g_json_indent_unit;
  auto first_valid_field = true;
  //First field is always END - ignore
  for(auto i=1u;i<query_config->get_num_queried_attributes();++i)
  {
    if(is_valid(i))
    {
      if(!first_valid_field)
        fptr << ",\n";
      fptr << indent_string << "\"" << (query_config->get_query_attribute_name(i)) << "\": ";
      if(i == ALT_query_idx)
        m_iterator->print_ALT(i, fptr);
      else
        m_iterator->print(i, fptr);
      first_valid_field = false;
    }
  }
  indent_string = indent_prefix+g_json_indent_unit;
  fptr << "\n" << indent_string << "}\n"<< indent_prefix << "}";
}

void GenomicsDBColumnarCell::print_csv(std::ostream& fptr, const VariantQueryConfig* query_config) const
{
  assert(query_config && query_config->is_bookkeeping_done());
  auto coords = get_coordinates();
  fptr << coords[0u] << "," << coords[1u];
  assert(query_config->is_defined_query_idx_for_known_field_enum(GVCF_END_IDX));
  auto end_query_idx = query_config->get_query_idx_for_known_field_enum(GVCF_END_IDX);
  auto ALT_query_idx = query_config->is_defined_query_idx_for_known_field_enum(GVCF_ALT_IDX)
    ? query_config->get_query_idx_for_known_field_enum(GVCF_ALT_IDX) : UINT64_MAX;
  auto end_position = *(reinterpret_cast<const int64_t*>(get_field_ptr_for_query_idx(end_query_idx)));
  fptr << "," << end_position;
  assert(coords[1] <= end_position);
  //First field is always END - ignore
  for(auto i=1ull;i<query_config->get_num_queried_attributes();++i)
  {
    fptr << ",";
    m_iterator->print_csv(i, fptr); //even if invalid, must print nulls
  }
  fptr << "\n";
}
