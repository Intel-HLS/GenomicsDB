#include "variant.h"

std::string g_non_reference_allele = "<NON_REF>";

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

//Variant functions
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
  //for(auto& ptr : m_fields)
  //ptr.reset(nullptr);
}

void Variant::print(std::ostream& fptr, const VariantQueryConfig* query_config) const
{
  fptr << "Interval:[ "<<m_col_begin <<", "<<m_col_end<<" ] Calls {";
  for(auto i=0ull;i<m_calls.size();++i)
  {
    fptr << " "<< i << " : {";
    m_calls[i].print(fptr, query_config ? query_config : m_query_config);
    fptr << " }";
  }
  fptr << " }\n";
}

void Variant::move_calls_to_separate_variants(std::vector<Variant>& variants, 
    std::vector<uint64_t>& query_row_idx_in_order)
{
  //Reverse order as gt_get_column uses reverse iterators
  for(int64_t i=query_row_idx_in_order.size()-1;i>=0;--i)
  {
    auto query_row_idx = query_row_idx_in_order[i];
    assert(query_row_idx < get_num_calls());
    variants.emplace_back(Variant());
    auto& curr_variant = variants[variants.size()-1u];
    auto& to_move_call = get_call(query_row_idx);
    //Set position of variant
    curr_variant.set_column_interval(to_move_call.get_column_begin(), to_move_call.get_column_end());
    curr_variant.add_call(std::move(to_move_call));
  }
}
