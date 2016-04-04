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

#include "variant.h"

auto json_indent_unit = "    ";

//GA4GHCallInfoToVariantIdx functions
bool GA4GHCallInfoToVariantIdx::find_or_insert(uint64_t begin, uint64_t end, const std::string& REF, 
        const std::vector<std::string>& ALT_vec, uint64_t& variant_idx)
{
  bool newly_inserted = true;
  auto col_begin_iter = m_begin_to_variant.find(begin);
  if(col_begin_iter == m_begin_to_variant.end())
    col_begin_iter = m_begin_to_variant.insert(std::pair<uint64_t, EndToVariantIdxTy>(begin, EndToVariantIdxTy())).first;
  auto& end_map = (*col_begin_iter).second;
  auto col_end_iter = end_map.find(end);
  if(col_end_iter == end_map.end())
    col_end_iter = end_map.insert(std::pair<uint64_t, REFToVariantIdxTy>(end, REFToVariantIdxTy())).first;
  auto& REF_map = (*col_end_iter).second;
  auto REF_iter = REF_map.find(REF);
  if(REF_iter == REF_map.end())
    REF_iter = REF_map.insert(std::pair<std::string, ALTSetToVariantIdxTy>(REF, ALTSetToVariantIdxTy())).first;
  auto& ALT_map = (*REF_iter).second;
  std::set<std::string> ALT_set;
  for(auto& ALT:ALT_vec)
    ALT_set.insert(ALT);
  auto ALT_iter = ALT_map.find(ALT_set);
  if(ALT_iter == ALT_map.end())
    ALT_iter = ALT_map.insert(std::pair<std::set<std::string>, uint64_t>(ALT_set, variant_idx)).first;
  else
  {
    newly_inserted = false;
    variant_idx = (*ALT_iter).second;
  }
  return newly_inserted;
}

bool GA4GHCallInfoToVariantIdx::find_or_insert(const VariantQueryConfig& query_config, VariantCall& to_move_call,
    uint64_t& variant_idx)
{
  const auto* REF_field_ptr =
    get_known_field_if_queried<VariantFieldString, true>(to_move_call, query_config, GVCF_REF_IDX); 
  const auto* ALT_field_ptr =
    get_known_field_if_queried<VariantFieldALTData, true>(to_move_call, query_config, GVCF_ALT_IDX);
  bool newly_inserted = true;
  //Checking for identical REF, ALT etc
  if(REF_field_ptr && ALT_field_ptr)
  {
    newly_inserted = find_or_insert(to_move_call.get_column_begin(), to_move_call.get_column_end(),
        REF_field_ptr->get(), ALT_field_ptr->get(), variant_idx);
  }
  return newly_inserted;
}

void GA4GHCallInfoToVariantIdx::clear()
{
  m_begin_to_variant.clear();
}

//GA4GHPagingInfo functions

void GA4GHPagingInfo::init_page_query()
{
  m_is_query_completed = false;
  m_last_row_idx = 0;
  m_last_column_idx = 0;
  m_num_handled_variants_in_last_column = 0u;
  m_num_variants_to_shift_left = 0u;
  m_num_variants_in_curr_page = 0u;
  deserialize_page_end();
}

void GA4GHPagingInfo::set_last_cell_info(std::vector<Variant>& variants,
    const uint64_t row_idx, const uint64_t column_idx, const unsigned num_last_column_variants_handled_after_curr_page)
{
#ifdef DEBUG
  assert(m_num_variants_in_curr_page == 0u);    //value set in init_page_query(), this function [resize] should never be called twice
  assert(m_num_handled_variants_in_last_column <= variants.size());
  assert(m_num_handled_variants_in_last_column + m_max_num_variants_per_page >= variants.size());
  assert(num_last_column_variants_handled_after_curr_page <= variants.size());
  //First set of variants must have same column as last_column_idx
  for(auto i=0u;i<m_num_handled_variants_in_last_column;++i)
    assert(variants[i].get_column_begin() == m_last_column_idx);
  //Last set of variants must have same column as column_idx
  for(auto i=variants.size()-num_last_column_variants_handled_after_curr_page;i<variants.size();++i)
    assert(variants[i].get_column_begin() == column_idx);
  //Curr page ends at same column as previous page, should have handled more variants [fwd progress]
  if(column_idx == m_last_column_idx)
  {
    assert(num_last_column_variants_handled_after_curr_page > m_num_handled_variants_in_last_column);
    assert(num_last_column_variants_handled_after_curr_page <= m_num_handled_variants_in_last_column
        + m_max_num_variants_per_page);
  }
#endif 
  m_num_variants_in_curr_page = variants.size() - m_num_handled_variants_in_last_column;
  m_last_row_idx = row_idx;
  m_last_column_idx = column_idx;
  m_num_variants_to_shift_left = m_num_handled_variants_in_last_column;
  m_num_handled_variants_in_last_column = num_last_column_variants_handled_after_curr_page;
}

void GA4GHPagingInfo::shift_left_variants(std::vector<Variant>& variants)
{
  //Remove variants handled in the previous page
  if(m_num_variants_to_shift_left)
  {
    //Shift left all variants in the vector
    auto num_variants_to_return = variants.size() - m_num_variants_to_shift_left;
    for(auto i=0u;i<num_variants_to_return;++i)
      variants[i] = std::move(variants[i+m_num_variants_to_shift_left]);
    variants.resize(num_variants_to_return);
  }
}

void GA4GHPagingInfo::serialize_page_end(const std::string& array_name)
{
  if(is_query_completed())
    m_last_page_end_token = "";
  else
    m_last_page_end_token = array_name + "_"
      + std::to_string(m_last_row_idx) + "_" 
      + std::to_string(m_last_column_idx) + "_"
      + std::to_string(m_num_handled_variants_in_last_column);
}

void GA4GHPagingInfo::deserialize_page_end()
{
  if(m_last_page_end_token == "")
  {
    m_last_column_idx = m_last_row_idx = 0ull;
    m_num_handled_variants_in_last_column = 0u;
    return;
  }
  char* dup_string = strdup(m_last_page_end_token.c_str());
  std::string row_string = "";
  std::string column_string = "";
  std::string num_handled_variants_string = "";
  char* saveptr = 0;
  char* ret_ptr = strtok_r(dup_string, "_", &saveptr);
  auto num_tokens = 0u;
  //Since array name may contain delimiters, tokenize and get the last 2 tokens only for row,col
  while(ret_ptr)
  {
    row_string = column_string;
    column_string = num_handled_variants_string;
    num_handled_variants_string = std::move(std::string(ret_ptr));
    ++num_tokens;
    ret_ptr = strtok_r(0, "_", &saveptr);
  }
  free(dup_string);
  if(num_tokens < 4u)   //<array_name>_<row>_<column>_<#variants_handled>
    throw InvalidGA4GHPageTokenException("Invalid GA4GH page token "+m_last_page_end_token+", TileDB-GA4GH page token should be of the form: <array_name>_<row>_<column>_<#handled_variants>");
  if(column_string.length() == 0u || row_string.length() == 0u || num_handled_variants_string.length() == 0u)
    throw InvalidGA4GHPageTokenException("Invalid GA4GH page token "+m_last_page_end_token+", TileDB-GA4GH page token should be of the form: <array_name>_<row>_<column>_<#handled_variants>");
  saveptr = 0;
  m_last_row_idx = strtoull(row_string.c_str(), &saveptr, 10);
  //Invalid number 
  if(saveptr == 0 ||  saveptr == row_string.c_str())
    throw InvalidGA4GHPageTokenException("Invalid GA4GH page token "+m_last_page_end_token+", TileDB-GA4GH page token should be of the form: <array_name>_<row>_<column> - row idx not detected");
  saveptr = 0;
  m_last_column_idx = strtoull(column_string.c_str(), &saveptr, 10);
  //Invalid number 
  if(saveptr == 0 ||  saveptr == column_string.c_str())
    throw InvalidGA4GHPageTokenException("Invalid GA4GH page token "+m_last_page_end_token+", TileDB-GA4GH page token should be of the form: <array_name>_<row>_<column> - column idx not detected");
  m_num_handled_variants_in_last_column = strtoull(num_handled_variants_string.c_str(), &saveptr, 10);
  //Invalid number 
  if(saveptr == 0 ||  saveptr == num_handled_variants_string.c_str())
    throw InvalidGA4GHPageTokenException("Invalid GA4GH page token "+m_last_page_end_token+", TileDB-GA4GH page token should be of the form: <array_name>_<row>_<column>_<#handled_variants> - #handled_variants not detected");
}


//Common functions used by VariantCall and Variant

//Function that copies field from src to dst
//Tries to minimize #re-allocations on the heap
void copy_field(std::unique_ptr<VariantFieldBase>& dst, const std::unique_ptr<VariantFieldBase>& src)
{
  auto dst_non_null = dst.get() ? 1u : 0u;
  auto src_non_null = src.get() ? 1u : 0u;
  auto combined = (dst_non_null << 1u) | src_non_null;
  switch(combined)
  {
    case 0u:     //both null, nothing to do
      break;
    case 1u:     //dst null, src non-null, create copy
      dst = std::move(std::unique_ptr<VariantFieldBase>(src->create_copy()));
      break;
    case 2u:     //dst non-null, src null, invalidate
      dst->set_valid(false);
      break;
    case 3u:     //both non-null, do fast copy
      dst->copy_from(src.get());
      break;
    default:
      break;
  }
}

void copy_fields(std::vector<std::unique_ptr<VariantFieldBase>>& dst, const std::vector<std::unique_ptr<VariantFieldBase>>& src)
{
  dst.resize(src.size());
  for(auto i=0u;i<src.size();++i)
    copy_field(dst[i], src[i]);
}

//VariantCall functions
void VariantCall::print(std::ostream& fptr, const VariantQueryConfig* query_config,
    const std::string& indent_prefix) const
{
  auto indent_string = indent_prefix+json_indent_unit;
  if(m_is_initialized && m_is_valid)
  {
    fptr << indent_prefix << "{\n";
    fptr << indent_string << "\"row\": "<<m_row_idx << ",\n";
    fptr << indent_string << "\"interval\": [ "<< m_col_begin << ", "<<m_col_end << " ],\n";
    fptr << indent_string << "\"fields\": {\n";
    indent_string += json_indent_unit;
    unsigned idx = 0u;
    auto first_valid_field = true;
    for(const auto& field : m_fields)
    {
      if(field.get() && field->is_valid())  //non null, valid field
      {
        if(!first_valid_field)
          fptr << ",\n";
        if(query_config)
          fptr << indent_string << "\"" << (query_config->get_query_attribute_name(idx)) << "\": ";
        else
          fptr << indent_string << "\"field_"<<idx<<"\": ";
        field->print(fptr);
        first_valid_field = false;
      }
      ++idx;
    }
    indent_string = indent_prefix+json_indent_unit;
    fptr << "\n" << indent_string << "}\n"<< indent_prefix << "}";
  }
}

void VariantCall::print_Cotton_JSON(std::ostream& fptr, unsigned field_idx) const
{
  if(m_is_initialized && m_is_valid)
  {
    assert(field_idx < m_fields.size());
    auto& field = m_fields[field_idx];
    if(field.get() && field->is_valid())  //non null, valid field
      field->print_Cotton_JSON(fptr);
    else
      fptr << "null";
  }
}

void VariantCall::reset_for_new_interval()
{
  m_is_initialized = false;
  m_is_valid = false;
  m_contains_deletion = false;
  //for(auto& ptr : m_fields)
  //ptr.reset(nullptr);
}

void VariantCall::copy_simple_members(const VariantCall& other)
{
  m_is_valid = other.is_valid();
  m_is_initialized = other.is_initialized();
  m_contains_deletion = other.m_contains_deletion;
  m_row_idx = other.get_row_idx();
  m_col_begin = other.m_col_begin;
  m_col_end = other.m_col_end;
}

/*
 * Performs move from other object
 */
void VariantCall::move_in(VariantCall& other)
{
  copy_simple_members(other);
  m_fields.resize(other.get_all_fields().size());
  unsigned idx = 0u;
  for(auto& other_field : other.get_all_fields())
  {
    set_field(idx, other_field);
    ++idx;
  }
}
/*
 * Creates copy of Call object
 */
void VariantCall::copy_from_call(const VariantCall& other)
{
  copy_simple_members(other);
  copy_fields(m_fields, other.get_all_fields());
}

void VariantCall::deep_copy_simple_members(const VariantCall& other)
{
  copy_simple_members(other);  //make copy
  resize(other.get_num_fields());
}

void VariantCall::binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const
{
  uint64_t add_size = 0ull;
  //is_valid, is_initialized, contains_deletion, row_idx, col_begin, col_end, num fields[unsigned]
  add_size = 3*sizeof(bool) + 3*sizeof(uint64_t) + sizeof(unsigned);
  RESIZE_BINARY_SERIALIZATION_BUFFER_IF_NEEDED(buffer, offset, add_size);
  //is_valid
  *(reinterpret_cast<bool*>(&(buffer[offset]))) = m_is_valid;
  offset += sizeof(bool);
  //is_initialized
  *(reinterpret_cast<bool*>(&(buffer[offset]))) = m_is_initialized;
  offset += sizeof(bool);
  //contains_deletion
  *(reinterpret_cast<bool*>(&(buffer[offset]))) = m_contains_deletion;
  offset += sizeof(bool);
  //row idx
  *(reinterpret_cast<uint64_t*>(&(buffer[offset]))) = m_row_idx;
  offset += sizeof(uint64_t);
  //column begin
  *(reinterpret_cast<uint64_t*>(&(buffer[offset]))) = m_col_begin;
  offset += sizeof(uint64_t);
  //column end
  *(reinterpret_cast<uint64_t*>(&(buffer[offset]))) = m_col_end;
  offset += sizeof(uint64_t);
  //num fields
  *(reinterpret_cast<unsigned*>(&(buffer[offset]))) = m_fields.size();
  offset += sizeof(unsigned);
  for(auto& field : m_fields)
  {
    //flag to represent valid field
    add_size = sizeof(bool);
    RESIZE_BINARY_SERIALIZATION_BUFFER_IF_NEEDED(buffer, offset, add_size);
    //is valid
    *(reinterpret_cast<bool*>(&(buffer[offset]))) = (field.get() && field->is_valid()) ? true : false;
    offset += sizeof(bool);
    //serialize valid field
    if(field.get() && field->is_valid())
      field->binary_serialize(buffer, offset);
  }
}

void VariantCall::binary_deserialize_header(const std::vector<uint8_t>& buffer, uint64_t& offset)
{
  //is_valid, is_initialized, contains_deletion, row_idx, col_begin, col_end, num fields[unsigned]
  assert(offset + 3*sizeof(bool) + 3*sizeof(uint64_t) + sizeof(unsigned) <= buffer.size());
  //is_valid
  m_is_valid = *(reinterpret_cast<const bool*>(&(buffer[offset])));
  offset += sizeof(bool);
  //is_initialized
  m_is_initialized = *(reinterpret_cast<const bool*>(&(buffer[offset])));
  offset += sizeof(bool);
  //contains_deletion
  m_contains_deletion = *(reinterpret_cast<const bool*>(&(buffer[offset])));
  offset += sizeof(bool);
  //row idx
  m_row_idx = *(reinterpret_cast<const uint64_t*>(&(buffer[offset])));
  offset += sizeof(uint64_t);
  //column begin
  m_col_begin = *(reinterpret_cast<const uint64_t*>(&(buffer[offset]))); 
  offset += sizeof(uint64_t);
  //column end
  m_col_end = *(reinterpret_cast<const uint64_t*>(&(buffer[offset])));
  offset += sizeof(uint64_t);
  //num fields
  auto num_fields = *(reinterpret_cast<const unsigned*>(&(buffer[offset])));
  offset += sizeof(unsigned);
  resize(num_fields);
}

//Variant functions
//FIXME: still assumes that Calls are allocated once and re-used across queries, need not be true
void Variant::reset_for_new_interval()
{
  for(auto& call : m_calls)
    call.reset_for_new_interval();
}

void Variant::resize_based_on_query()
{
  assert(m_query_config);
  assert(m_query_config->is_bookkeeping_done());
  //Initialize VariantCall vector and pointer vector
  uint64_t num_rows = m_query_config->get_num_rows_to_query();
  resize(num_rows, m_query_config->get_num_queried_attributes());
  for(uint64_t i=0ull;i<num_rows;++i)
  {
    uint64_t row_idx = m_query_config->get_array_row_idx_for_query_row_idx(i);
    m_calls[i].set_row_idx(row_idx);
  }
}

void Variant::print(std::ostream& fptr, const VariantQueryConfig* query_config, const std::string& indent_prefix) const
{
  fptr << indent_prefix << "{\n";
  std::string indent_string = indent_prefix+json_indent_unit;
  fptr << indent_string << "\"interval\": [ "<<m_col_begin <<", "<<m_col_end<<" ],\n";
  fptr << indent_string << " \"common_fields\" : {\n";
  indent_string += json_indent_unit;
  auto idx = 0u;
  auto first_valid_field = true;
  for(const auto& field : m_fields)
  {
    if(field.get() && field->is_valid())  //non null, valid field
    {
      if(!first_valid_field)
        fptr << ",\n";
      fptr << indent_string;
      if(query_config)
        fptr << "\"" << (query_config->get_query_attribute_name(m_common_fields_query_idxs[idx])) << "\": ";
      else
        fptr << "\"field_" << idx << "\": ";
      field->print(fptr);
      first_valid_field = false;
    }
    ++idx;
  }
  indent_string = indent_prefix + json_indent_unit;
  fptr << "\n" << indent_string <<"},\n";
  fptr << indent_string << "\"variant_calls\": [\n";
  indent_string += json_indent_unit;
  auto call_idx = 0ull;
  for(auto iter=begin();iter!=end();++iter)
  {
    if(call_idx > 0ull)
      fptr << ",\n";
    (*iter).print(fptr, query_config ? query_config : m_query_config, indent_string);
    ++call_idx;
  }
  indent_string = indent_prefix + json_indent_unit;
  fptr << "\n" << indent_string << "]\n";
  fptr << indent_prefix << "}";
}

void print_field(std::ostream& fptr,
                const std::vector<Variant>& variants,
                const std::vector<uint64_t>& starts,
                const std::vector<uint64_t>& ends,
                const ContigInfo* contig_info,
                const unsigned& field_type,
                unsigned attribute_index = 0) 
{
  auto first_valid = true;
  int64_t tile_position;

  for (size_t worker_idx = 0; worker_idx < starts.size(); ++worker_idx)
  {
    for (uint64_t i = starts[worker_idx]; i < ends[worker_idx]; ++i)
    {
      const auto& variant = variants[i];
      for(const auto& call : variant)
      {
        if(!first_valid)
          fptr << ",";
        else
          first_valid = false;
        switch(field_type) {
          case INDICES_IDX:
            fptr << call.get_row_idx();
            break;
          case START_IDX:
            tile_position = call.get_column_begin();
            assert(tile_position >= contig_info->m_tiledb_column_offset);
            assert((tile_position - contig_info->m_tiledb_column_offset) < contig_info->m_length);
            // Adding 1 as TileDB is 0-based and genomics (VCF) is 1-based
            fptr << (tile_position - contig_info->m_tiledb_column_offset + 1);
            break;
          case END_IDX:
            tile_position = call.get_column_end();
            assert(tile_position >= contig_info->m_tiledb_column_offset);
            assert((tile_position - contig_info->m_tiledb_column_offset) < contig_info->m_length);
            // Adding 1 as TileDB is 0-based and genomics (VCF) is 1-based
            fptr << (tile_position - contig_info->m_tiledb_column_offset + 1);
            break;
          case ATTRIBUTE_IDX:
            call.print_Cotton_JSON(fptr, attribute_index);
            break;
          default:
            assert(0);
        }
      }
    }
  }
}

void print_fields(std::ostream& fptr,
                  const std::vector<Variant>& variants,
                  const VariantQueryConfig& query_config,
                  const std::vector<uint64_t>& starts,
                  const std::vector<uint64_t>& ends,
                  const ContigInfo* contig_info) 
{
  std::string indent = json_indent_unit;
  fptr << indent + "\"indices\" : [ ";
  print_field(fptr, variants, starts, ends, contig_info, INDICES_IDX);
  fptr << " ],\n";

  fptr << indent + "\"POSITION\" : [ ";
  print_field(fptr, variants, starts, ends, contig_info, START_IDX);
  fptr << " ],\n";

  fptr << indent + "\"END\" : [ ";
  print_field(fptr, variants, starts, ends, contig_info, END_IDX);
  fptr << " ],\n";

  //other attributes, start from 1 as the first queried attribute is always END
  for(auto i=1u;i<query_config.get_num_queried_attributes();++i)
  {
    fptr << indent + "\"" + query_config.get_query_attribute_name(i) + "\" : [ ";
    print_field(fptr, variants, starts, ends, contig_info, ATTRIBUTE_IDX, i);
    fptr << " ]";
    if(i+1u >= query_config.get_num_queried_attributes())       //last query, no comma
      fptr << "\n";
    else
      fptr << ",\n";
  }
}

void print_Cotton_JSON(std::ostream& fptr, const std::vector<Variant>& variants, const VariantQueryConfig& query_config, const VidMapper* id_mapper)
{
  assert(query_config.is_bookkeeping_done());
  std::string contig_name;
  int64_t contig_position;
  ContigInfo contig_info;

  // Currently Cotton-JSON is called only on single position queries
  if (!id_mapper->get_contig_location(query_config.get_column_begin(0), contig_name, contig_position))
    throw VidMapperException("print_and_update_contig_position: Invalid position ");
  // Initialize contig_info
  if (!id_mapper->get_contig_info(contig_name, contig_info))
    throw VidMapperException("print_and_update_contig_position: Invalid contig name : " + contig_name );

  fptr << "{\n";
  std::vector<uint64_t> starts(1, 0);
  std::vector<uint64_t> ends(1, variants.size());
  print_fields(fptr, variants, query_config, starts, ends, &contig_info);
  fptr << "}\n";
}

void print_and_update_contig_position(std::ostream& fptr,
                                      const int64_t& query_start,
                                      const int64_t& query_end,
                                      int64_t& start_position,
                                      int64_t& end_position,
                                      std::string& contig_name,
                                      std::string& prev_contig_name,
                                      ContigInfo* contig_info,
                                      const VidMapper* id_mapper)
{
  if (!(id_mapper->get_contig_location(query_start, contig_name, start_position) &&
           id_mapper->get_contig_location(query_end, contig_name, end_position) ))
  {
    throw VidMapperException("print_and_update_contig_position: Invalid position ");
  }
  // Adding 1 as TileDB is 0-based and genomics (VCF) is 1-based
  ++start_position;
  ++end_position;

  // Check if the contig has changed then close the braces and add the contig
  if (prev_contig_name.size() == 0)
  {
    // This is the first time we will be printing contig information
    fptr << "\"" << contig_name << "\": {\n";
    prev_contig_name = contig_name;
    // Initialize contig_info
    if (!id_mapper->get_contig_info(contig_name, *contig_info))
    {
      throw VidMapperException("print_and_update_contig_position: Invalid contig name : " + contig_name );
    }
  }
  else if (prev_contig_name.compare(contig_name) != 0)
  {
    fptr << "},\n";
    fptr << "\"" << contig_name << "\": {\n";
    prev_contig_name = contig_name;
    // Update contig_info for the new contig
    if (!id_mapper->get_contig_info(contig_name, *contig_info))
    {
      throw VidMapperException("print_and_update_contig_position: Invalid contig name : " + contig_name );
    }
  }
  else {
    // We are in the same contig as before
    // Add a , and output the next queried position in query
    fptr << ",\n";
  }
}

void print_positions_json_split_by_column(std::ostream& fptr,
                                          const std::vector<Variant>& variants,
                                          const std::vector<uint64_t>& query_column_lengths,
                                          const std::vector<uint64_t>& num_column_intervals,
                                          const std::vector<uint64_t>& queried_column_positions,
                                          const VariantQueryConfig& query_config,
                                          const VidMapper* id_mapper)
{
  unsigned queried_index = 0u;
  unsigned variant_start_index = 0u;
  unsigned variant_start_offset = 0u;
  ContigInfo contig_info;
  std::string contig_name;
  std::string prev_contig_name;
  int64_t start_position;
  int64_t end_position;
  std::vector<uint64_t> ends(1, 0ull);
  std::vector<uint64_t> starts(1, 0ull);
  // Start the JSON
  fptr << "{\n";
  for (auto interval = 0ull; interval < num_column_intervals.size(); ++interval)
  {
    for (auto i = 0ull; i < num_column_intervals[interval]; ++i, ++queried_index)
    {
      print_and_update_contig_position(fptr,
                                      queried_column_positions[queried_index * 2],
                                      queried_column_positions[queried_index * 2 + 1],
                                      start_position, end_position,
                                      contig_name, prev_contig_name,
                                      &contig_info, id_mapper
                                      );
      starts[0] = variant_start_index;
      ends[0] = variant_start_offset + query_column_lengths[queried_index];
      fptr << "\"" << start_position;
      if (start_position != end_position)
      {
        fptr  << "_" << end_position;
      }
      fptr << "\" : {\n";
      print_fields(fptr, variants, query_config, starts, ends, &contig_info);
      // Close data for current query position
      fptr << "}\n";
      // set start index to the end variant index
      variant_start_index = ends[0];
      // If this is the last queried column from the current rank/worker,
      // then update offset to the last but one printed variant index
      if (i == (num_column_intervals[interval] - 1))
        variant_start_offset += query_column_lengths[queried_index];
    }
  }
  // Close the Contig level
  fptr << "}\n";
  // Close the JSON
  fptr << "}\n";
}

void print_positions_json_split_by_row(std::ostream& fptr,
                                      const std::vector<Variant>& variants,
                                      const VariantQueryConfig& query_config,
                                      const std::vector<uint64_t>& query_column_lengths,
                                      const std::vector<uint64_t>& num_column_intervals,
                                      const VidMapper* id_mapper) {
  // This function works only when the tile db instances are split by samples
  assert(query_config.is_bookkeeping_done());
  auto num_partitions = num_column_intervals.size();
  auto num_queried_columns = query_config.get_num_column_intervals();
  std::vector<uint64_t> ends(num_partitions, 0ull);
  std::vector<uint64_t> starts(num_partitions, 0ull);
  ContigInfo contig_info;
  std::string contig_name;
  std::string prev_contig_name;
  int64_t start_position;
  int64_t end_position;
  fptr << "{\n";
  for (auto query_column_idx = 0ull; query_column_idx < num_queried_columns; ++query_column_idx)
  {
    auto query_start = query_config.get_column_begin(query_column_idx);
    auto query_end   = query_config.get_column_end(query_column_idx);
    print_and_update_contig_position(fptr, query_start, query_end, start_position, end_position,
                              contig_name, prev_contig_name, &contig_info, id_mapper);
    uint64_t start  = 0;
    uint64_t aggregated_worker_offset = 0;
    for (auto worker_idx = 0ull; worker_idx < num_partitions; ++worker_idx)
    {
      // length assumes start is at 0, but with multiple workers start is offset,
      // so add the total length of variants from the previous worker
      ends[worker_idx] = start + query_column_lengths[(worker_idx * num_queried_columns) + query_column_idx];
      // Add Offset to start idx with the length from previous column
      if (query_column_idx > 0)
      {
        start += query_column_lengths[(worker_idx * num_queried_columns) + (query_column_idx - 1)];
      }
      starts[worker_idx] = start;
      // Update start to the start of the variant idx for the next worker
      aggregated_worker_offset += query_column_lengths[num_column_intervals[worker_idx] - 1];
      start = aggregated_worker_offset;
    }
    fptr << "\"" << start_position;
    if (start_position != end_position)
    {
      fptr  << "_" << end_position;
    }
    fptr << "\" : {\n";
    print_fields(fptr, variants, query_config, starts, ends, &contig_info);
    // Close data for current query position
    fptr << "}\n";
  }
  // Close the Contig level
  fptr << "}\n";
  // Close the JSON
  fptr << "}\n";
}

void Variant::get_column_sorted_call_idx_vec(std::vector<uint64_t>& query_row_idx_in_order)
{
  //no ordering exists in query_row_idx_in_order - do a column major sort first
  query_row_idx_in_order.resize(get_num_calls());
  //iterate over valid calls
  auto valid_call_idx = 0ull;
  for(auto iter=begin();iter!=end();++iter,++valid_call_idx)
    query_row_idx_in_order[valid_call_idx] = iter.get_call_idx_in_variant();
  query_row_idx_in_order.resize(valid_call_idx);
  std::sort(query_row_idx_in_order.begin(), query_row_idx_in_order.end(), VariantCallIdxColumnMajorLT(this));
}

void Variant::move_calls_to_separate_variants(const VariantQueryConfig& query_config, std::vector<Variant>& variants, 
    std::vector<uint64_t>& query_row_idx_in_order, GA4GHCallInfoToVariantIdx& call_info_2_variant, GA4GHPagingInfo* paging_info)
{
#ifdef DUPLICATE_CELL_AT_END
  get_column_sorted_call_idx_vec(query_row_idx_in_order);
#else
  if(query_row_idx_in_order.size() == 0u)
    return;
#endif
  uint64_t last_column_idx = paging_info ? paging_info->get_last_column() : 0u;
  auto num_last_column_variants_handled_after_curr_page = 0u;
  bool stop_inserting_new_variants = false;
#ifdef DUPLICATE_CELL_AT_END
  //Sorted in column major order - so normal order
  for(auto i=0ull;i<query_row_idx_in_order.size();++i)
#else
  //Reverse order as gt_get_column uses reverse iterators
  for(int64_t i=query_row_idx_in_order.size()-1;i>=0;--i)
#endif
  {
    auto query_row_idx = query_row_idx_in_order[i];
    assert(query_row_idx < get_num_calls());
    auto& to_move_call = get_call(query_row_idx);
    auto curr_row_idx = to_move_call.get_row_idx();
    auto curr_column_idx = to_move_call.get_column_begin();
    //If this is a continued query, only return results after the last page
    if(paging_info && paging_info->handled_previously(curr_row_idx, curr_column_idx))
      continue;
    auto newly_inserted = move_call_to_variant_vector(query_config, to_move_call, variants, call_info_2_variant,
        stop_inserting_new_variants); 
    //If paging, complex logic for checking page end 
    PAGE_END_CHECK_LOGIC 
    last_column_idx = curr_column_idx;
  }
}

void Variant::copy_simple_members(const Variant& other)
{
  m_query_config = other.m_query_config;
  m_col_begin = other.m_col_begin;
  m_col_end = other.m_col_end;
}

//Memory de-allocation
void Variant::clear()
{
  for(auto& call : m_calls)
    call.clear();
  m_calls.clear();
  m_fields.clear();
  m_common_fields_query_idxs.clear();
}

//Function that moves information from other to self
void Variant::move_in(Variant& other)
{
  //Copy simple primitives
  copy_simple_members(other);
  //Move Calls
  m_calls.resize(other.get_num_calls());
  for(auto i=0ull;i<other.get_num_calls();++i)
    m_calls[i] = std::move(other.get_call(i));
  //Move common fields
  resize_common_fields(other.get_num_common_fields());
  for(auto i=0u;i<other.get_num_common_fields();++i)
    set_common_field(i, other.get_query_idx_for_common_field(i), other.get_common_field(i));
}

void Variant::copy_from_variant(const Variant& other)
{
  //Copy simple primitive members
  copy_simple_members(other);
  //Copy Calls
  m_calls.resize(other.get_num_calls());
  for(auto i=0ull;i<other.get_num_calls();++i)
    m_calls[i].copy_from_call(other.get_call(i));  //make copy
  copy_common_fields(other);
}

void Variant::copy_common_fields(const Variant& other)
{
  //Resize common fields
  resize_common_fields(other.get_num_common_fields());
  if(get_num_common_fields())
  {
    //Copy query idxs
    memcpy(&(m_common_fields_query_idxs[0]), &(other.get_common_fields_query_idxs()[0]), m_common_fields_query_idxs.size()*sizeof(unsigned));
    copy_fields(m_fields, other.get_common_fields());
  }
}

void Variant::deep_copy_simple_members(const Variant& other)
{
  //Copy simple primitive members
  copy_simple_members(other);
  //Copy simple primitive members of member calls
  m_calls.resize(other.get_num_calls());
  for(auto i=0ull;i<other.get_num_calls();++i)
    m_calls[i].deep_copy_simple_members(other.get_call(i));
}

void Variant::binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const
{
  uint64_t add_size = 0ull;
  //Header - column begin, end, num_calls, num_common_fields[unsigned]
  add_size = 3*sizeof(uint64_t) + sizeof(unsigned);
  RESIZE_BINARY_SERIALIZATION_BUFFER_IF_NEEDED(buffer, offset, add_size);
  //Col begin
  *(reinterpret_cast<uint64_t*>(&(buffer[offset]))) = m_col_begin;
  offset += sizeof(uint64_t);
  //Col end
  *(reinterpret_cast<uint64_t*>(&(buffer[offset]))) = m_col_end;
  offset += sizeof(uint64_t);
  //num calls
  *(reinterpret_cast<uint64_t*>(&(buffer[offset]))) = get_num_calls();
  offset += sizeof(uint64_t);
  //num common fields
  *(reinterpret_cast<unsigned*>(&(buffer[offset]))) = get_num_common_fields();
  offset += sizeof(unsigned);
  //Serialize calls
  for(auto i=0ull;i<get_num_calls();++i)
    m_calls[i].binary_serialize(buffer, offset);
  //Common fields
  for(auto i=0u;i<get_num_common_fields();++i)
  {
    //Flag representing whether common field is valid and query idx for common field
    add_size = sizeof(bool) + sizeof(unsigned);
    RESIZE_BINARY_SERIALIZATION_BUFFER_IF_NEEDED(buffer, offset, add_size);
    auto& curr_field = m_fields[i];
    //Flag representing whether common field is not null
    *(reinterpret_cast<bool*>(&(buffer[offset]))) = (curr_field.get() && curr_field->is_valid()) ? true : false;
    offset += sizeof(bool);
    //Query idx for common field
    *(reinterpret_cast<unsigned*>(&(buffer[offset]))) = get_query_idx_for_common_field(i);
    offset += sizeof(unsigned);
    //Serialize common field
    if(curr_field.get() && curr_field->is_valid())
      curr_field->binary_serialize(buffer, offset);
  }
}

void Variant::binary_deserialize_header(const std::vector<uint8_t>& buffer, uint64_t& offset, unsigned num_queried_attributes)
{
  assert(offset + 3*sizeof(uint64_t) + sizeof(unsigned) <= buffer.size());
  //Col begin
  auto col_begin = *(reinterpret_cast<const uint64_t*>(&(buffer[offset])));
  offset += sizeof(uint64_t);
  //Col end
  auto col_end = *(reinterpret_cast<const uint64_t*>(&(buffer[offset])));
  offset += sizeof(uint64_t);
  //num calls
  auto num_calls = *(reinterpret_cast<const uint64_t*>(&(buffer[offset])));
  offset += sizeof(uint64_t);
  //num common fields
  auto num_common_fields = *(reinterpret_cast<const unsigned*>(&(buffer[offset])));
  offset += sizeof(unsigned);
  //column info
  set_column_interval(col_begin, col_end);
  //Resize variant
  resize(num_calls, num_queried_attributes);
  resize_common_fields(num_common_fields);
}

bool move_call_to_variant_vector(const VariantQueryConfig& query_config, VariantCall& to_move_call,
    std::vector<Variant>& variants, GA4GHCallInfoToVariantIdx& call_info_2_variant, bool stop_inserting_new_variants)
{
  uint64_t variant_idx = variants.size();
  bool newly_inserted = call_info_2_variant.find_or_insert(query_config, to_move_call, variant_idx);
  if(newly_inserted && !stop_inserting_new_variants)
    variants.emplace_back(Variant());
  //variant_idx can be >= variants.size() if stop_inserting_new_variants is true
  assert(variant_idx < variants.size() || stop_inserting_new_variants);
  if((!newly_inserted && variant_idx < variants.size()) || !stop_inserting_new_variants)
  {
    Variant& curr_variant = variants[variant_idx];
    //Set position of variant
    curr_variant.set_column_interval(to_move_call.get_column_begin(), to_move_call.get_column_end());
    curr_variant.add_call(std::move(to_move_call));
  }
  return newly_inserted;
}

void print_variants(const std::vector<Variant>& variants,
                    const std::string& output_format,
                    const VariantQueryConfig& query_config,
                    std::ostream& fptr,
                    const bool is_partitioned_by_column,
                    const VidMapper* id_mapper,
                    const std::vector<uint64_t>& query_column_lengths,
                    const std::vector<uint64_t>& num_column_intervals,
                    const std::vector<uint64_t>& queried_column_positions,
                    bool output_directly)
{
  static const std::unordered_map<std::string, unsigned> format_2_enum = {
    { "Cotton-JSON", COTTON_JSON_OUTPUT_FORMAT_IDX },
    { "Positions-JSON", POSITIONS_JSON_OUTPUT_FORMAT_IDX },
    { "GA4GH", GA4GH_OUTPUT_FORMAT_IDX }
  };
  unsigned output_format_idx = DEFAULT_OUTPUT_FORMAT_IDX;
  auto iter = format_2_enum.find(output_format);
  if(iter != format_2_enum.end())
    output_format_idx = (*iter).second;
  std::stringstream ss;
  //if output stream is a stringstream, use it directly, else output to stringstream first
  std::ostream& optr = (output_directly || dynamic_cast<std::ostringstream*>(&fptr)) ? fptr : ss;
  optr << std::fixed;
  optr << std::setprecision(6);
  switch(output_format_idx)
  {
    case COTTON_JSON_OUTPUT_FORMAT_IDX:
      print_Cotton_JSON(optr, variants, query_config, id_mapper);
      break;
    case POSITIONS_JSON_OUTPUT_FORMAT_IDX:
      if (is_partitioned_by_column)
        print_positions_json_split_by_column(optr, variants,
                                            query_column_lengths,
                                            num_column_intervals,
                                            queried_column_positions,
                                            query_config,
                                            id_mapper);
      else
        print_positions_json_split_by_row(optr, variants, query_config, query_column_lengths, num_column_intervals, id_mapper);
      break;
    case DEFAULT_OUTPUT_FORMAT_IDX:
    default:
      {
        optr << "{\n" << json_indent_unit << "\"variants\": [\n";
        auto variant_idx = 0ull;
        auto indent_prefix = std::string(json_indent_unit)+json_indent_unit;
        for(const auto& variant : variants)
        {
          if(variant_idx > 0ull)
            optr << ",\n";
          variant.print(optr, &query_config, indent_prefix);
          ++variant_idx;
        }
        optr << "\n"<< json_indent_unit << "]\n";
        optr << "}\n";
        break;
      }
  }
  //If using stringstream as a temp buffer, print out to fptr
  if(&optr == &ss)
  {
    std::string buffer;
#define STRINGSTREAM_BUFFER_SIZE 65536u
    buffer.resize(STRINGSTREAM_BUFFER_SIZE);
    while(!(ss.eof()) && !(ss.fail()))
    {
      ss.read(&(buffer[0]), STRINGSTREAM_BUFFER_SIZE);
      auto count = ss.gcount();
      fptr.write(&(buffer[0]), count);
    }
    ss.clear();
  }
}
