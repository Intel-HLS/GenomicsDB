#include "variant.h"

std::string g_non_reference_allele = "<NON_REF>";

GTColumn::GTColumn(int64_t col, uint64_t row_num) {
  ALT_.resize(row_num);
  col_ = col;
  REF_.resize(row_num);
  PL_.resize(row_num);
  AF_.resize(row_num);
  AN_.resize(row_num);
  AC_.resize(row_num);
}

void GTColumn::reset()
{
  for(auto i=0ull;i<REF_.size();++i)
  {
    REF_[i].clear();
    ALT_[i].clear();
    PL_[i].clear();
    AF_.clear();
    AN_.clear();
    AC_.clear();
  }
}
//FIXME: still assumes that Calls are allocated once and re-used across queries, need not be true
void Variant::reset_for_new_interval()
{
  for(auto& call : m_calls)
    call.reset_for_new_interval();
}

void VariantCall::reset_for_new_interval()
{
  m_is_initialized = false;
  m_is_valid = false;
  for(auto& ptr : m_fields)
    ptr.reset(nullptr);
}

void Variant::print(std::ostream& fptr) const
{
  fptr << "Interval:[ "<<m_col_begin <<", "<<m_col_end<<" ] Calls {";
  for(auto i=0ull;i<m_calls.size();++i)
  {
    fptr << " "<< i << " : {";
    m_calls[i].print(fptr, m_query_config);
    fptr << " }";
  }
  fptr << " }\n";
}

void VariantCall::print(std::ostream& fptr, const VariantQueryConfig* query_config) const
{
  if(m_is_initialized && m_is_valid)
  {
    fptr << " row : "<<m_row_idx << ", ";
    unsigned idx = 0u;
    for(const auto& field : m_fields)
    {
      if(field.get())  //non null field
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
