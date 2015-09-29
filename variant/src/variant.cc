#include "variant.h"

std::string g_non_reference_allele = "<NON_REF>";
uint32_t bcf_float_missing    = 0x7F800001;
uint32_t bcf_float_vector_end = 0x7F800002;

fi_pair bcf_float_missing_union = { .i = bcf_float_missing };

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
    {
      field->print_Cotton_JSON(fptr);
      return;
    }
    else
      fptr << "null";
  }
}

void VariantCall::reset_for_new_interval()
{
  m_is_initialized = false;
  m_is_valid = false;
  //for(auto& ptr : m_fields)
  //ptr.reset(nullptr);
}

void VariantCall::copy_simple_members(const VariantCall& other)
{
  m_is_valid = other.is_valid();
  m_is_initialized = other.is_initialized();
  m_row_idx = other.get_row_idx();
  m_col_begin = other.m_col_begin;
  m_col_end = other.m_col_end;
}

/*
 * Performs move from other object
 */
void VariantCall::move_in(VariantCall& other)
{
  clear();
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
  clear();
  copy_simple_members(other);
  m_fields.resize(other.get_all_fields().size());
  unsigned idx = 0u;
  for(const auto& other_field : other.get_all_fields())
  {
    set_field(idx, other_field.get() ? other_field->create_copy() : 0); //if non-null, create copy, else null
    ++idx;
  }
}

void VariantCall::binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const
{
  uint64_t add_size = 0ull;
  //is_valid, is_initialized, row_idx, col_begin, col_end, num fields
  add_size = 2*sizeof(bool) + 3*sizeof(uint64_t) + sizeof(unsigned);
  if(offset + add_size > buffer.size())
    buffer.resize(offset + add_size + 1024ull); //extra space
  //is_valid
  *(reinterpret_cast<bool*>(&(buffer[offset]))) = m_is_valid;
  offset += sizeof(bool);
  //is_initialized
  *(reinterpret_cast<bool*>(&(buffer[offset]))) = m_is_initialized;
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
    //flag to represent non-null field
    add_size = sizeof(uint8_t);
    if(offset + add_size > buffer.size())
      buffer.resize(offset + add_size + 1024ull); //extra space
    //is non-null
    buffer[offset] = field.get() ? 1u : 0u;
    offset += sizeof(uint8_t);
    //serialize field
    if(field.get())
      field->binary_serialize(buffer, offset);
  }
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
    if(field.get())  //non null field
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
  bool first = true;
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
    std::vector<uint64_t>& query_row_idx_in_order, GA4GHCallInfoToVariantIdx& call_info_2_variant)
{
  if(query_row_idx_in_order.size() == 0)
    return;

  //Reverse order as gt_get_column uses reverse iterators
  for(int64_t i=query_row_idx_in_order.size()-1;i>=0;--i)
  {
    auto query_row_idx = query_row_idx_in_order[i];
    assert(query_row_idx < get_num_calls());
    auto& to_move_call = get_call(query_row_idx);
    move_call_to_variant_vector(query_config, to_move_call, variants, call_info_2_variant);
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
  //De-allocates existing data
  clear();
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
  //De-allocates existing data
  clear();
  //Copy simple primitive members
  copy_simple_members(other);
  //Copy Calls
  m_calls.resize(other.get_num_calls());
  for(auto i=0ull;i<other.get_num_calls();++i)
    m_calls[i].copy_from_call(other.get_call(i));  //make copy
  //Copy common fields
  resize_common_fields(other.get_num_common_fields());
  for(auto i=0u;i<other.get_num_common_fields();++i)
    set_common_field(i, other.get_query_idx_for_common_field(i), 
        other.get_common_field(i).get() ? other.get_common_field(i)->create_copy() : 0);    //copy if non-null, else null
}

void Variant::binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const
{
  uint64_t add_size = 0ull;
  //Header - column begin, end, num_calls, num_common_fields[unsigned]
  add_size = 3*sizeof(uint64_t) + sizeof(unsigned);
  if(offset + add_size > buffer.size())
    buffer.resize(offset + add_size + 1024u);   //large size
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
  for(auto i=0u;i<get_num_common_fields();++i)
  {
    //Flag representing whether common field is not null and query idx for common field
    add_size = sizeof(uint8_t) + sizeof(unsigned);
    if(offset + add_size > buffer.size())
      buffer.resize(offset+1024u);      //larger size
    //Flag representing whether common field is not null
    buffer[offset] = m_fields[i].get() ? 1 : 0;
    offset += sizeof(uint8_t);
    //Query idx for common field
    *(reinterpret_cast<unsigned*>(&(buffer[offset]))) = get_query_idx_for_common_field(i);
    offset += sizeof(unsigned);
    //Serialize common field
    if(m_fields[i].get())
      m_fields[i]->binary_serialize(buffer, offset);
  }
}

void move_call_to_variant_vector(const VariantQueryConfig& query_config, VariantCall& to_move_call,
    std::vector<Variant>& variants, GA4GHCallInfoToVariantIdx& call_info_2_variant)
{
  uint64_t variant_idx = variants.size();
  bool newly_inserted = call_info_2_variant.find_or_insert(query_config, to_move_call, variant_idx);
  if(newly_inserted)
    variants.emplace_back(Variant());
  Variant& curr_variant = variants[variant_idx];
  //Set position of variant
  curr_variant.set_column_interval(to_move_call.get_column_begin(), to_move_call.get_column_end());
  curr_variant.add_call(std::move(to_move_call));
}
