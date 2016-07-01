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

#include "variant_operations.h"

template<>
std::string get_zero_value() { return ""; }

/*
  Copied from defunct gamgee library
  Remaps data dependent on number of genotypes to the new order of alleles as specified in alleles_LUT
  @param input_data - vector of data values for a given input sample as stored in TileDB
  @param input_call_idx 
  @param alleles_LUT LUT mapping alleles from input to merged alleles list
  @param num_merged_alleles
  @param NON_REF_exists
  @param alt_alleles_only flag that determines whether only the ALT alleles will be used or all alleles
  @param remapped_data - data structure in which remapped info will be stored 
  @num_calls_with_valid_data - keeps track of how many samples had valid values for given genotype idx
 */
template<class DataType>
void  VariantOperations::remap_data_based_on_alleles(const std::vector<DataType>& input_data,
    const uint64_t input_call_idx, 
    const CombineAllelesLUT& alleles_LUT, const unsigned num_merged_alleles, bool NON_REF_exists, bool alt_alleles_only,
    RemappedDataWrapperBase& remapped_data,
    std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value) {
  //index of NON_REF in merged variant
  const auto merged_non_reference_allele_idx = NON_REF_exists ? 
    static_cast<int64_t>(static_cast<int>(num_merged_alleles-1)) : lut_missing_value;
  //index of NON_REF in input sample
  const auto input_non_reference_allele_idx = NON_REF_exists ? 
    alleles_LUT.get_input_idx_for_merged(input_call_idx, merged_non_reference_allele_idx) : lut_missing_value;
  //Loop over alleles - only ALT or all alleles (BCF_VL_A or BCF_VL_R)
  unsigned length = alt_alleles_only ? num_merged_alleles-1u: num_merged_alleles;
  for (auto j=0u;j<length;++j) {
    auto allele_j = alt_alleles_only ?  j+1u : j;
    auto input_j_allele = alleles_LUT.get_input_idx_for_merged(input_call_idx, allele_j);
    if (CombineAllelesLUT::is_missing_value(input_j_allele))	//no mapping found for current allele in input gvcf
    {
      if(CombineAllelesLUT::is_missing_value(input_non_reference_allele_idx))	//input did not have NON_REF allele
      {
        *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, j))) = (missing_value);
        continue;
      }
      else //input contains NON_REF allele, use its idx
        input_j_allele = input_non_reference_allele_idx;
    }
    assert(!alt_alleles_only || input_j_allele > 0u);   //if only ALT alleles are used, then input_j_allele must be non-0
    auto input_j = alt_alleles_only ? input_j_allele-1u : input_j_allele;
    *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, j))) = 
      input_data[input_j];
    ++(num_calls_with_valid_data[j]);
  }
}
/*
  Copied from defunct gamgee library
  Remaps data dependent on number of genotypes to the new order of alleles as specified in alleles_LUT
  @param input_data - vector of data values for a given input sample as stored in TileDB
  @param input_call_idx 
  @param alleles_LUT LUT mapping alleles from input to merged alleles list
  @param num_merged_alleles  
  @param remapped_data - data structure in which remapped info will be stored 
  @num_calls_with_valid_data - keeps track of how many samples had valid values for given genotype idx
 */
template<class DataType>
void  VariantOperations::remap_data_based_on_genotype(const std::vector<DataType>& input_data,
    const uint64_t input_call_idx, 
    const CombineAllelesLUT& alleles_LUT, const unsigned num_merged_alleles, bool NON_REF_exists,
    RemappedDataWrapperBase& remapped_data,
    std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value) {
  //index of NON_REF in merged variant
  const auto merged_non_reference_allele_idx = NON_REF_exists ? 
    static_cast<int64_t>(static_cast<int>(num_merged_alleles-1)) : lut_missing_value;
  //index of NON_REF in input sample
  const auto input_non_reference_allele_idx = NON_REF_exists ? 
    alleles_LUT.get_input_idx_for_merged(input_call_idx, merged_non_reference_allele_idx) : lut_missing_value;
  //Loop over all possible genotype combinations
  for (auto allele_j = 0u; allele_j < num_merged_alleles; ++allele_j) {
    auto input_j_allele = alleles_LUT.get_input_idx_for_merged(input_call_idx, allele_j);
    if (CombineAllelesLUT::is_missing_value(input_j_allele))	//no mapping found for current allele in input gvcf
    {
      if(CombineAllelesLUT::is_missing_value(input_non_reference_allele_idx))	//input did not have NON_REF allele
      {
	//fill in missing values for all genotypes with allele_j as one component
	for(auto allele_k = allele_j; allele_k < num_merged_alleles;++allele_k)
        {
          auto gt_idx = bcf_alleles2gt(allele_j, allele_k);
	  *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, gt_idx))) = (missing_value);
        }
	continue;	//skip to next value of allele_j
      }
      else //input contains NON_REF allele, use its idx
	input_j_allele = input_non_reference_allele_idx;
    }
    for (auto allele_k = allele_j; allele_k < num_merged_alleles; ++allele_k) {
      auto gt_idx = bcf_alleles2gt(allele_j, allele_k);
      auto input_k_allele = alleles_LUT.get_input_idx_for_merged(input_call_idx, allele_k);
      if (CombineAllelesLUT::is_missing_value(input_k_allele))	//no mapping found for current allele in input gvcf
      {
	if(CombineAllelesLUT::is_missing_value(input_non_reference_allele_idx))	//input did not have NON_REF allele
	{
	  *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, gt_idx))) = (missing_value);
	  continue;	//skip to next value of allele_k
	}
	else //input has NON_REF, use its idx
	  input_k_allele = input_non_reference_allele_idx;
      }
      *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, gt_idx))) = 
        input_data[bcf_alleles2gt(input_j_allele, input_k_allele)];
      ++(num_calls_with_valid_data[gt_idx]);
    }
  }
}

//Variant handler functions
template<class DataType>
void VariantFieldHandler<DataType>::remap_vector_data(std::unique_ptr<VariantFieldBase>& orig_field_ptr, uint64_t curr_call_idx_in_variant, 
    const CombineAllelesLUT& alleles_LUT, unsigned num_merged_alleles, bool non_ref_exists,
    unsigned length_descriptor, unsigned num_merged_elements, RemappedVariant& remapper_variant)
{
  auto* raw_orig_field_ptr = orig_field_ptr.get();
  if(raw_orig_field_ptr == 0)
    return;
  //Assert that ptr is of type VariantFieldPrimitiveVectorData<DataType>
  assert(dynamic_cast<VariantFieldPrimitiveVectorData<DataType>*>(raw_orig_field_ptr));
  auto* orig_vector_field_ptr = static_cast<VariantFieldPrimitiveVectorData<DataType>*>(raw_orig_field_ptr);
  //Resize and zero vector
  m_num_calls_with_valid_data.resize(num_merged_elements);
  memset(&(m_num_calls_with_valid_data[0]), 0, num_merged_elements*sizeof(uint64_t));
  /*Remap field in copy (through remapper_variant)*/
  if(KnownFieldInfo::is_length_descriptor_genotype_dependent(length_descriptor))
    VariantOperations::remap_data_based_on_genotype<DataType>( 
        orig_vector_field_ptr->get(), curr_call_idx_in_variant, 
        alleles_LUT, num_merged_alleles, non_ref_exists, 
        remapper_variant, m_num_calls_with_valid_data, m_bcf_missing_value); 
  else 
    VariantOperations::remap_data_based_on_alleles<DataType>( 
        orig_vector_field_ptr->get(), curr_call_idx_in_variant, 
        alleles_LUT, num_merged_alleles, non_ref_exists, 
        KnownFieldInfo::is_length_descriptor_only_ALT_alleles_dependent(length_descriptor),
        remapper_variant, m_num_calls_with_valid_data, m_bcf_missing_value);
}

template<class DataType>
bool VariantFieldHandler<DataType>::get_valid_median(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, void* output_ptr)
{
  m_median_compute_vector.resize(variant.get_num_calls());
  auto valid_idx = 0u;
  //Iterate over valid calls
  for(auto iter=variant.begin(), end_iter = variant.end();iter != end_iter;++iter)
  {
    auto& curr_call = *iter;
    auto& field_ptr = curr_call.get_field(query_idx);
    //Valid field
    if(field_ptr.get() && field_ptr->is_valid())
    {
      //Must always be vector<DataType>
      auto* ptr = dynamic_cast<VariantFieldPrimitiveVectorData<DataType>*>(field_ptr.get());
      assert(ptr); 
      assert((ptr->get()).size() > 0u);
      auto val = ptr->get()[0u];
      if(is_bcf_valid_value<DataType>(val))
        m_median_compute_vector[valid_idx++] = val;
    }
  }
  if(valid_idx == 0u)   //no valid fields found
    return false;
  auto mid_point = valid_idx/2u;
  std::nth_element(m_median_compute_vector.begin(), m_median_compute_vector.begin()+mid_point, m_median_compute_vector.begin()+valid_idx);
  auto result_ptr = reinterpret_cast<DataType*>(output_ptr);
  *result_ptr = m_median_compute_vector[mid_point];
  return true;
}

template<class DataType>
bool VariantFieldHandler<DataType>::get_valid_sum(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, void* output_ptr)
{
  DataType sum = get_zero_value<DataType>();
  auto valid_idx = 0u;
  //Iterate over valid calls
  for(auto iter=variant.begin(), end_iter = variant.end();iter != end_iter;++iter)
  {
    auto& curr_call = *iter;
    auto& field_ptr = curr_call.get_field(query_idx);
    //Valid field
    if(field_ptr.get() && field_ptr->is_valid())
    {
      //Must always be vector<DataType>
      auto* ptr = dynamic_cast<VariantFieldPrimitiveVectorData<DataType>*>(field_ptr.get());
      assert(ptr); 
      assert((ptr->get()).size() > 0u);
      auto val = ptr->get()[0u];
      if(is_bcf_valid_value<DataType>(val))
      {
        sum += val;
        ++valid_idx;
      }
    }
  }
  if(valid_idx == 0u)   //no valid fields found
    return false;
  auto result_ptr = reinterpret_cast<DataType*>(output_ptr);
  *result_ptr = sum;
  return true;
}

template<class DataType>
bool VariantFieldHandler<DataType>::compute_valid_element_wise_sum(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, void* output_ptr, unsigned num_elements)
{
  DataType* sum = reinterpret_cast<DataType*>(output_ptr);
  for(auto i=0u;i<num_elements;++i)
    sum[i] = get_zero_value<DataType>();
  auto valid_idx = 0u;
  //Iterate over valid calls
  for(auto iter=variant.begin(), end_iter = variant.end();iter != end_iter;++iter)
  {
    auto& curr_call = *iter;
    auto& field_ptr = curr_call.get_field(query_idx);
    //Valid field
    if(field_ptr.get() && field_ptr->is_valid())
    {
      //Must always be vector<DataType>
      auto* ptr = dynamic_cast<VariantFieldPrimitiveVectorData<DataType>*>(field_ptr.get());
      assert(ptr);
      auto& vec = ptr->get();
      assert(vec.size() > 0u);
      for(auto i=0ull;i<vec.size() && i<num_elements;++i)
      {
        auto val = vec[i];
        if(is_bcf_valid_value<DataType>(val))
        {
          sum[i] += val;
          ++valid_idx;
        }
      }
    }
  }
  return (valid_idx > 0u);
}

template<class DataType>
bool VariantFieldHandler<DataType>::collect_and_extend_fields(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, const void ** output_ptr, unsigned& num_elements,
        const bool use_missing_values_only_not_vector_end, const bool use_vector_end_only)
{
  auto max_elements_per_call = 0u;
  auto valid_idx = 0u;
  //Iterate over valid calls and obtain max fields over all calls
  for(auto iter=variant.begin(), end_iter = variant.end();iter != end_iter;++iter)
  {
    auto& curr_call = *iter;
    auto& field_ptr = curr_call.get_field(query_idx);
    //Valid field
    if(field_ptr.get() && field_ptr->is_valid())
    {
      max_elements_per_call = std::max<unsigned>(max_elements_per_call, field_ptr->length());
      ++valid_idx;
    }
  }
  if(valid_idx == 0u)   //no valid fields found
    return false;
  //Resize the extended field vector
  if(variant.get_num_calls()*max_elements_per_call > m_extended_field_vector.size())
    m_extended_field_vector.resize(variant.get_num_calls()*max_elements_per_call);
  auto extended_field_vector_idx = 0u;
  //Iterate over all calls, invalid calls also
  for(auto call_idx=0ull;call_idx<variant.get_num_calls();++call_idx)
  {
    auto& curr_call = variant.get_call(call_idx);
    auto& field_ptr = curr_call.get_field(query_idx);
    //#elements inserted for this call
    auto num_elements_inserted = 0u;
    //Valid field in a valid call
    if(curr_call.is_valid() && field_ptr.get() && field_ptr->is_valid())
    {
      assert(field_ptr->get_raw_pointer());
      memcpy(&(m_extended_field_vector[extended_field_vector_idx]), field_ptr->get_raw_pointer(), field_ptr->length()*sizeof(DataType));
      num_elements_inserted = field_ptr->length();
      extended_field_vector_idx += num_elements_inserted;
    }
    if(num_elements_inserted == 0u) //no elements inserted, insert missing value first
    {
      m_extended_field_vector[extended_field_vector_idx] = use_vector_end_only ? get_bcf_vector_end_value<DataType>()
        : get_bcf_missing_value<DataType>();
      ++num_elements_inserted;
      ++extended_field_vector_idx;
    }
    //Pad with vector end values, handles invalid fields also
    //Except when producing records for htsjdk BCF2 - htsjdk has no support for vector end values
    auto padded_value = use_missing_values_only_not_vector_end ? get_bcf_missing_value<DataType>()
        : get_bcf_vector_end_value<DataType>();
    for(;num_elements_inserted<max_elements_per_call;++num_elements_inserted,++extended_field_vector_idx)
      m_extended_field_vector[extended_field_vector_idx] = padded_value;
  }
  assert(extended_field_vector_idx <= m_extended_field_vector.size());
  *output_ptr = reinterpret_cast<const void*>(&(m_extended_field_vector[0]));
  num_elements = extended_field_vector_idx;
  return true;
}
//Explicit template instantiation
template class VariantFieldHandler<int>;
template class VariantFieldHandler<unsigned>;
template class VariantFieldHandler<int64_t>;
template class VariantFieldHandler<uint64_t>;
template class VariantFieldHandler<float>;
template class VariantFieldHandler<double>;
template class VariantFieldHandler<std::string>;
template class VariantFieldHandler<char>;
