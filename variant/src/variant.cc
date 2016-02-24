#include "variant.h"

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
void VariantCall::print(std::ostream& fptr, const VariantQueryConfig* query_config) const
{
  if(m_is_initialized && m_is_valid)
  {
    fptr << " row : "<<m_row_idx << ", ";
    fptr << "interval : [ "<< m_col_begin << ", "<<m_col_end << " ], ";
    unsigned idx = 0u;
    for(const auto& field : m_fields)
    {
      if(field.get() && field->is_valid())  //non null, valid field
      {
        if(query_config)
          fptr << (query_config->get_query_attribute_name(idx)) << " : ";
        field->print(fptr);
        fptr << ", ";
      }
      ++idx;
    }
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

void Variant::print(std::ostream& fptr, const VariantQueryConfig* query_config) const
{
  fptr << "Interval:[ "<<m_col_begin <<", "<<m_col_end<<" ]";
  fptr << " Common fields : { ";
  auto idx = 0u;
  for(const auto& field : m_fields)
  {
    if(field.get() && field->is_valid())  //non null, valid field
    {
      if(query_config)
        fptr << (query_config->get_query_attribute_name(m_common_fields_query_idxs[idx])) << " : ";
      field->print(fptr);
      fptr << ", ";
    }
    ++idx;
  }
  fptr <<" } Calls {";
  for(auto i=0ull;i<m_calls.size();++i)
  {
    fptr << " "<< i << " : {";
    m_calls[i].print(fptr, query_config ? query_config : m_query_config);
    fptr << " }";
  }
  fptr << " }\n";
}

void print_Cotton_JSON(std::ostream& fptr, const std::vector<Variant>& variants, const VariantQueryConfig& query_config)
{
  assert(query_config.is_bookkeeping_done());
  fptr << "{\n";
  std::string indent = "    ";
  //Row idxs
  {
    auto first_valid = true;
    fptr << indent + "\"indices\" : [ ";
    for(const auto& variant : variants)
    {
      for(const auto& call : variant)
      {
        if(!first_valid)
          fptr << "," << call.get_row_idx();
        else
        {
          fptr << call.get_row_idx();
          first_valid = false;
        }
      }
    }
    fptr << " ],\n";
  }
  //Start
  {
    auto first_valid = true;
    fptr << indent + "\"start\" : [ ";
    for(const auto& variant : variants)
    {
      for(const auto& call : variant)
      {
        if(!first_valid)
          fptr << "," << call.get_column_begin();
        else
        {
          fptr << call.get_column_begin();
          first_valid = false;
        }
      }
    }
    fptr << " ],\n";
  }
  //END
  {
    auto first_valid = true;
    fptr << indent + "\"end\" : [ ";
    for(const auto& variant : variants)
    {
      for(const auto& call : variant)
      {
        if(!first_valid)
          fptr << "," << call.get_column_end();
        else
        {
          fptr << call.get_column_end();
          first_valid = false;
        }
      }
    }
    fptr << " ],\n";
  }
  //other attributes, start from 1 as the first queried attribute is always END
  for(auto i=1u;i<query_config.get_num_queried_attributes();++i)
  {
    auto first_valid = true;
    fptr << indent + "\"" + query_config.get_query_attribute_name(i) + "\" : [ ";
    for(const auto& variant : variants)
    {
      for(const auto& call : variant)
      {
        if(!first_valid)
        {
          fptr << ",";
          call.print_Cotton_JSON(fptr, i);
        }
        else
        {
          call.print_Cotton_JSON(fptr, i);
          first_valid = false;
        }
      }
    }
    fptr << " ]";
    if(i+1u >= query_config.get_num_queried_attributes())       //last query, no comma
      fptr << "\n";
    else
      fptr << ",\n";
  }
  fptr << "}\n";
}

void Variant::move_calls_to_separate_variants(const VariantQueryConfig& query_config, std::vector<Variant>& variants, 
    std::vector<uint64_t>& query_row_idx_in_order, GA4GHCallInfoToVariantIdx& call_info_2_variant, GA4GHPagingInfo* paging_info)
{
  if(query_row_idx_in_order.size() == 0u)
    return;
  uint64_t last_column_idx = paging_info ? paging_info->get_last_column() : 0u;
  auto num_last_column_variants_handled_after_curr_page = 0u;
  bool stop_inserting_new_variants = false;
  //Reverse order as gt_get_column uses reverse iterators
  for(int64_t i=query_row_idx_in_order.size()-1;i>=0;--i)
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

void print_variants(const std::vector<Variant>& variants, const std::string& output_format, const VariantQueryConfig& query_config,
    std::ostream& fptr, bool output_directly)
{
  static const std::unordered_map<std::string, unsigned> format_2_enum = {
    { "Cotton-JSON", COTTON_JSON_OUTPUT_FORMAT_IDX },
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
      print_Cotton_JSON(optr, variants, query_config);
      break;
    case DEFAULT_OUTPUT_FORMAT_IDX:
    default:
      for(const auto& variant : variants)
        variant.print(optr, &query_config);
      break;
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
