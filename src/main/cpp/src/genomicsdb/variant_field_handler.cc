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
    //Input data could have been truncated due to missing values - if so, put missing value
    if(static_cast<size_t>(input_j) >= input_data.size())
      *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, j))) = missing_value;
    else
    {
      *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, j))) =
        input_data[input_j];
      if(is_bcf_valid_value<DataType>(input_data[input_j]))
        ++(num_calls_with_valid_data[j]);
    }
  }
}

template<class DataType>
void  VariantOperations::remap_data_based_on_genotype_haploid(const std::vector<DataType>& input_data,
    const uint64_t input_call_idx,
    const CombineAllelesLUT& alleles_LUT,
    const unsigned num_merged_alleles, bool NON_REF_exists,
    RemappedDataWrapperBase& remapped_data,
    std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value) {
  //index of NON_REF in merged variant
  const auto merged_non_reference_allele_idx = NON_REF_exists ? 
    static_cast<int64_t>(static_cast<int>(num_merged_alleles-1)) : lut_missing_value;
  //index of NON_REF in input sample
  const auto input_non_reference_allele_idx = NON_REF_exists ? 
    alleles_LUT.get_input_idx_for_merged(input_call_idx, merged_non_reference_allele_idx) : lut_missing_value;
  //Loop over all possible genotype combinations == #alleles
  for (auto allele_j = 0u; allele_j < num_merged_alleles; ++allele_j) {
    auto input_j_allele = alleles_LUT.get_input_idx_for_merged(input_call_idx, allele_j);
    auto gt_idx = allele_j;
    if (CombineAllelesLUT::is_missing_value(input_j_allele))	//no mapping found for current allele in input gvcf
    {
      if(CombineAllelesLUT::is_missing_value(input_non_reference_allele_idx))	//input did not have NON_REF allele
      {
        *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, gt_idx))) = (missing_value);
        continue;	//skip to next value of allele_j
      }
      else //input contains NON_REF allele, use its idx
        input_j_allele = input_non_reference_allele_idx;
    }
    //Input data could have been truncated due to missing values - if so, put missing value
    if(static_cast<size_t>(input_j_allele) >= input_data.size())
      *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, gt_idx))) = (missing_value);
    else
    {
      auto input_gt_idx = input_j_allele;
      *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, gt_idx))) =
        input_data[input_gt_idx];
      if(is_bcf_valid_value<DataType>(input_data[input_gt_idx]))
        ++(num_calls_with_valid_data[gt_idx]);
    }
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
void  VariantOperations::remap_data_based_on_genotype_diploid(const std::vector<DataType>& input_data,
    const uint64_t input_call_idx,
    const CombineAllelesLUT& alleles_LUT,
    const unsigned num_merged_alleles, bool NON_REF_exists,
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
      auto input_gt_idx = (bcf_alleles2gt(input_j_allele, input_k_allele));
      //Input data could have been truncated due to missing values - if so, put missing value
      if(static_cast<size_t>(input_gt_idx) >= input_data.size())
        *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, gt_idx))) = (missing_value);
      else
      {
        *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, gt_idx))) =
          input_data[input_gt_idx];
        if(is_bcf_valid_value<DataType>(input_data[input_gt_idx]))
          ++(num_calls_with_valid_data[gt_idx]);
      }
    }
  }
}


/*
 * Reorders fields whose length and order depend on the number of genotypes (BCF_VL_G)
 * for general ploidy
 */
template<class DataType>
void VariantOperations::remap_data_based_on_genotype_general(const std::vector<DataType>& input_data,
    const uint64_t input_call_idx,
    const CombineAllelesLUT& alleles_LUT,
    const unsigned num_merged_alleles, bool NON_REF_exists, const unsigned ploidy,
    RemappedDataWrapperBase& remapped_data,
    std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value,
    std::vector<int>& remapped_allele_idx_vec_for_current_gt_combination,
    std::vector<std::pair<int, int> >& ploidy_index_allele_index_stack,
    std::vector<int>& input_call_allele_idx_vec_for_current_gt_combination,
    remap_operator_function_type<DataType> op
)
{
  if(ploidy == 0u)
    return;
  //index of NON_REF in merged variant
  const auto merged_non_reference_allele_idx = NON_REF_exists ? 
    static_cast<int64_t>(static_cast<int>(num_merged_alleles-1)) : lut_missing_value;
  //index of NON_REF in input sample
  const auto input_non_reference_allele_idx = NON_REF_exists ? 
    alleles_LUT.get_input_idx_for_merged(input_call_idx, merged_non_reference_allele_idx) : lut_missing_value;
#define SET_PLOIDY_INDEX_IN_STACK_ELEMENT(X,Y) ((X).first = (Y))
#define SET_ALLELE_INDEX_IN_STACK_ELEMENT(X,Y) ((X).second = (Y))
#define GET_PLOIDY_INDEX_IN_STACK_ELEMENT(X) ((X).first)
#define GET_ALLELE_INDEX_IN_STACK_ELEMENT(X) ((X).second)
  remapped_allele_idx_vec_for_current_gt_combination.resize(ploidy+1u); //+1 to avoid unnecessary if statements in the while loop
  input_call_allele_idx_vec_for_current_gt_combination.resize(ploidy);
  //Enumerate genotypes based on method described in 
  //http://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
  //Use custom "stack" instead of a recursive function call
  //"top" of the stack is the last element of the vector
  //resize to max #genotypes - avoids frequent dynamic memory allocations/frees
  ploidy_index_allele_index_stack.resize(KnownFieldInfo::get_number_of_genotypes(num_merged_alleles-1u, ploidy));
  //In each iteration, generate all genotypes where the last ploidy corresponds to the allele
  //corresponding to top_level_allele_idx
  auto allele_idx = 0;
  auto ploidy_idx = 0;
  //In each iteration of the loop, let the top of the stack contain (P=x, A=y)
  //The subsequent iterations will generate all genotype combinations corresponding to
  //gt[x] == y first before moving on to elements in the stack below
  //Initializer element in the stack set (P=ploidy, A=#alleles-1)
  //Thus, the while loop will generate all genotype combinations for ploidies 0..ploidy-1
  //with alleles 0..alleles-1
  SET_PLOIDY_INDEX_IN_STACK_ELEMENT(ploidy_index_allele_index_stack[0u], ploidy);
  SET_ALLELE_INDEX_IN_STACK_ELEMENT(ploidy_index_allele_index_stack[0u], num_merged_alleles-1);
  auto num_elements_in_stack = 1u;
  auto remapped_gt_idx = 0ull;
  while(num_elements_in_stack > 0u)
  {
    auto& top_stack_element = ploidy_index_allele_index_stack[num_elements_in_stack-1u];
    allele_idx = GET_ALLELE_INDEX_IN_STACK_ELEMENT(top_stack_element);
    ploidy_idx = GET_PLOIDY_INDEX_IN_STACK_ELEMENT(top_stack_element);
    --num_elements_in_stack; //popped stack
    assert(ploidy_idx >= 0 && static_cast<size_t>(ploidy_idx) < remapped_allele_idx_vec_for_current_gt_combination.size());
    remapped_allele_idx_vec_for_current_gt_combination[ploidy_idx] = allele_idx;
    //Assigned one allele idx for all ploidys
    if(ploidy_idx == 0)
    {
      auto curr_genotype_combination_contains_missing_allele_for_input = false;
      for(auto i=0u;i<ploidy;++i)
      {
        auto input_allele_idx = alleles_LUT.get_input_idx_for_merged(input_call_idx, remapped_allele_idx_vec_for_current_gt_combination[i]);
        if (CombineAllelesLUT::is_missing_value(input_allele_idx))	//no mapping found for current allele in input gvcf
        {
          input_call_allele_idx_vec_for_current_gt_combination[i] = input_non_reference_allele_idx; //set to NON_REF idx, possibly missing
          curr_genotype_combination_contains_missing_allele_for_input =
            curr_genotype_combination_contains_missing_allele_for_input ||
            CombineAllelesLUT::is_missing_value(input_non_reference_allele_idx);
        }
        else
          input_call_allele_idx_vec_for_current_gt_combination[i] = input_allele_idx;
      }
      op(input_data,
          input_call_idx,
          alleles_LUT,
          num_merged_alleles, NON_REF_exists,
          curr_genotype_combination_contains_missing_allele_for_input,
          ploidy,
          remapped_data,
          num_calls_with_valid_data, missing_value,
          remapped_allele_idx_vec_for_current_gt_combination,
          remapped_gt_idx,
          input_call_allele_idx_vec_for_current_gt_combination);
      ++remapped_gt_idx;
    }
    else
    {
      --ploidy_idx; //current ploidy_idx
      //Reverse order so that alleles with lower idx are closer to the top of the stack
      for(auto i=allele_idx;i>=0;--i)
      {
        assert(num_elements_in_stack < ploidy_index_allele_index_stack.size());
        auto& curr_stack_element = ploidy_index_allele_index_stack[num_elements_in_stack];
        SET_PLOIDY_INDEX_IN_STACK_ELEMENT(curr_stack_element, ploidy_idx);
        SET_ALLELE_INDEX_IN_STACK_ELEMENT(curr_stack_element, i);
        ++num_elements_in_stack;
      }
    }
  }
}

uint64_t VariantOperations::get_genotype_index(std::vector<int>& allele_idx_vec, const bool is_sorted)
{
  //To get genotype combination index for the input, alleles must be in sorted order
  if(!is_sorted)
    std::sort(allele_idx_vec.begin(), allele_idx_vec.end());
  switch(allele_idx_vec.size()) //ploidy
  {
    case 0u:
      return 0u;
    case 1u:
      return allele_idx_vec[0u];
    case 2u:
      return bcf_alleles2gt(allele_idx_vec[0u], allele_idx_vec[1u]);
    default:
      {
        auto gt_idx = 0ull;
        //From http://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
        for(auto i=0ull;i<allele_idx_vec.size();++i)
          gt_idx += VariantOperations::nCr(i+allele_idx_vec[i], allele_idx_vec[i]-1);
        return gt_idx;
      }
  }
}

template<class DataType>
void VariantOperations::reorder_field_based_on_genotype_index(const std::vector<DataType>& input_data,
    const uint64_t input_call_idx,
    const CombineAllelesLUT& alleles_LUT,
    const unsigned num_merged_alleles, bool NON_REF_exists,
    const bool curr_genotype_combination_contains_missing_allele_for_input,
    const unsigned ploidy,
    RemappedDataWrapperBase& remapped_data,
    std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value,
    const std::vector<int>& remapped_allele_idx_vec_for_current_gt_combination,
    const uint64_t remapped_gt_idx,
    std::vector<int>& input_call_allele_idx_vec_for_current_gt_combination
    )
{
  if(curr_genotype_combination_contains_missing_allele_for_input) //no genotype in input corresponding to this allele combination
  {
    //Put missing value
    *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, remapped_gt_idx))) = (missing_value);
    return;
  }
  auto input_gt_idx = VariantOperations::get_genotype_index(input_call_allele_idx_vec_for_current_gt_combination, false);
  assert(remapped_gt_idx < KnownFieldInfo::get_number_of_genotypes(num_merged_alleles-1u, ploidy));
  //Input data could have been truncated due to missing values - if so, put missing value
  if(input_gt_idx >= input_data.size())
    *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, remapped_gt_idx))) = (missing_value);
  else
  {
    *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, remapped_gt_idx))) =
      input_data[input_gt_idx];
    if(is_bcf_valid_value<DataType>(input_data[input_gt_idx]))
      ++(num_calls_with_valid_data[remapped_gt_idx]);
  }
}

//Wrapper function
template<class DataType>
void  VariantOperations::remap_data_based_on_genotype(const std::vector<DataType>& input_data,
    const uint64_t input_call_idx,
    const CombineAllelesLUT& alleles_LUT,
    const unsigned num_merged_alleles, bool NON_REF_exists, const unsigned ploidy,
    RemappedDataWrapperBase& remapped_data,
    std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value,
    std::vector<int>& remapped_allele_idx_vec_for_current_gt_combination, std::vector<std::pair<int, int> >& ploidy_index_allele_index_stack,
    std::vector<int>& input_call_allele_idx_vec_for_current_gt_combination)
{
  switch(ploidy)
  {
    case 1u:
      remap_data_based_on_genotype_haploid(input_data,
          input_call_idx,
          alleles_LUT,
          num_merged_alleles, NON_REF_exists,
          remapped_data,
          num_calls_with_valid_data, missing_value);
      break;
    case 2u:
      remap_data_based_on_genotype_diploid(input_data,
          input_call_idx,
          alleles_LUT,
          num_merged_alleles, NON_REF_exists,
          remapped_data,
          num_calls_with_valid_data, missing_value);
      break;
    default:  //why not let this case handle diploid and haploid? A. speed, diploid is common case
      remap_data_based_on_genotype_general<DataType>(input_data,
          input_call_idx,
          alleles_LUT,
          num_merged_alleles, NON_REF_exists, ploidy,
          remapped_data,
          num_calls_with_valid_data, missing_value,
          remapped_allele_idx_vec_for_current_gt_combination, ploidy_index_allele_index_stack,
          input_call_allele_idx_vec_for_current_gt_combination,
          VariantOperations::reorder_field_based_on_genotype_index<DataType>);
      break;
  }
}

//Variant handler functions
template<class DataType>
void VariantFieldHandler<DataType>::remap_vector_data(std::unique_ptr<VariantFieldBase>& orig_field_ptr, uint64_t curr_call_idx_in_variant, 
    const CombineAllelesLUT& alleles_LUT,
    unsigned num_merged_alleles, bool non_ref_exists, const unsigned ploidy,
    const FieldLengthDescriptor& length_descriptor, unsigned num_merged_elements, RemappedVariant& remapper_variant)
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
  if(length_descriptor.is_length_genotype_dependent())
    VariantOperations::remap_data_based_on_genotype<DataType>( 
        orig_vector_field_ptr->get(), curr_call_idx_in_variant, 
        alleles_LUT,
        num_merged_alleles, non_ref_exists, ploidy,
        remapper_variant, m_num_calls_with_valid_data, m_bcf_missing_value,
        m_allele_idx_vec_for_current_genotype, m_ploidy_index_alleles_index_stack,
        m_input_call_allele_idx_vec); 
  else 
    VariantOperations::remap_data_based_on_alleles<DataType>( 
        orig_vector_field_ptr->get(), curr_call_idx_in_variant, 
        alleles_LUT, num_merged_alleles, non_ref_exists, 
        length_descriptor.is_length_only_ALT_alleles_dependent(),
        remapper_variant, m_num_calls_with_valid_data, m_bcf_missing_value);
}

template<class DataType>
bool VariantFieldHandler<DataType>::get_valid_median(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, void* output_ptr, unsigned& num_valid_elements)
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
        unsigned query_idx, void* output_ptr, unsigned& num_valid_elements)
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
  num_valid_elements = valid_idx;
  if(valid_idx == 0u)   //no valid fields found
    return false;
  auto result_ptr = reinterpret_cast<DataType*>(output_ptr);
  *result_ptr = sum;
  return true;
}


template<class DataType>
bool VariantFieldHandler<DataType>::get_valid_mean(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, void* output_ptr, unsigned& num_valid_elements)
{
  auto status = get_valid_sum(variant, query_config, query_idx, output_ptr, num_valid_elements);
  if(status)
  {
    auto result_ptr = reinterpret_cast<DataType*>(output_ptr);
    *result_ptr = (*result_ptr)/num_valid_elements;
  }
  return status;
}

//Mean is undefined for strings
template<>
bool VariantFieldHandler<std::string>::get_valid_mean(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, void* output_ptr, unsigned& num_valid_elements)
{
  throw VariantOperationException("Mean is an undefined operation for combining string fields");
  return false;
}

template<class DataType>
bool VariantFieldHandler<DataType>::compute_valid_element_wise_sum(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, const void** output_ptr, unsigned& num_elements)
{
  auto num_valid_elements = 0u;
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
      if(vec.size() > m_element_wise_operations_result.size())
        m_element_wise_operations_result.resize(vec.size());
      for(auto i=0ull;i<vec.size();++i)
      {
        auto val = vec[i];
        if(is_bcf_valid_value<DataType>(val))
        {
          if(i < num_valid_elements && is_bcf_valid_value<DataType>(m_element_wise_operations_result[i]))
            m_element_wise_operations_result[i] += val;
          else
          {
            m_element_wise_operations_result[i] = val;
            if(i >= num_valid_elements)
            {
              //Set all elements after the last valid value upto i to missing
              for(auto j=num_valid_elements;j<i;++j)
                m_element_wise_operations_result[j] = get_bcf_missing_value<DataType>();
              num_valid_elements = i+1u;
            }
          }
        }
      }
    }
  }
  if(num_valid_elements > 0u)
    m_element_wise_operations_result.resize(num_valid_elements);
  (*output_ptr) = &(m_element_wise_operations_result[0]);
  num_elements = num_valid_elements;
  return (num_valid_elements > 0u);
}

template<class DataType>
bool VariantFieldHandler<DataType>::concatenate_field(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, const void** output_ptr, unsigned& num_elements)
{
  auto curr_result_size = 0ull;
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
      if(curr_result_size + vec.size() > m_element_wise_operations_result.size())
        m_element_wise_operations_result.resize(curr_result_size+vec.size());
      memcpy(&(m_element_wise_operations_result[curr_result_size]), &(vec[0]), vec.size()*sizeof(DataType));
      curr_result_size += vec.size();
    }
  }
  if(curr_result_size > 0u)
    m_element_wise_operations_result.resize(curr_result_size);
  (*output_ptr) = &(m_element_wise_operations_result[0]);
  num_elements = curr_result_size;
  return (curr_result_size > 0u);
}

//I apologize - really should not have distinguished between vector<non-char> and string [vector<char>] :(
//In the columnar field branch, I hope to remove this unnecessary distinction
template<>
bool VariantFieldHandler<char>::concatenate_field(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, const void** output_ptr, unsigned& num_elements)
{
  auto curr_result_size = 0ull;
  //Iterate over valid calls
  for(auto iter=variant.begin(), end_iter = variant.end();iter != end_iter;++iter)
  {
    auto& curr_call = *iter;
    auto& field_ptr = curr_call.get_field(query_idx);
    //Valid field
    if(field_ptr.get() && field_ptr->is_valid())
    {
      auto* ptr = dynamic_cast<VariantFieldString*>(field_ptr.get());
      assert(ptr);
      auto& vec = ptr->get();
      if(curr_result_size + vec.size() > m_element_wise_operations_result.size())
        m_element_wise_operations_result.resize(curr_result_size+vec.size());
      memcpy(&(m_element_wise_operations_result[curr_result_size]), &(vec[0]), vec.size()*sizeof(char));
      curr_result_size += vec.size();
    }
  }
  if(curr_result_size > 0u)
    m_element_wise_operations_result.resize(curr_result_size);
  (*output_ptr) = &(m_element_wise_operations_result[0]);
  num_elements = curr_result_size;
  return (curr_result_size > 0u);
}

template<class DataType>
bool VariantFieldHandler<DataType>::collect_and_extend_fields(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, const void ** output_ptr, uint64_t& num_elements,
        const bool use_missing_values_only_not_vector_end, const bool use_vector_end_only,
        const bool is_GT_field)
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
    //WARNING: bunch of horrible hacks to deal with VCF spec and its implementations: htslib and htsjdk
    if(num_elements_inserted == 0u) //no elements inserted for this call, insert missing value first
    {
      //use_vector_end_only - true only for string fields when the Java interface is used
      //use_missing_values_only_not_vector_end - true only when the Java interface is used
      m_extended_field_vector[extended_field_vector_idx] =
        (use_vector_end_only || (is_GT_field && !use_missing_values_only_not_vector_end))
        ? get_bcf_vector_end_value<DataType>()
        : ((is_GT_field && use_missing_values_only_not_vector_end)
            ? get_bcf_gt_missing_value<DataType>() //why? htsjdk does not handle a record where GT is missing in 1 sample, present in another
            : get_bcf_missing_value<DataType>());
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

//The following function is called in variant_operations.cc also. For some compilers, explicit instantiation
//of VariantFieldHandler is insufficient. The following function must also be instantiated, else weird link time
//errors are seen for variant_operations.cc. I'm guessing there is some inline error, but cannot be sure.
template
void VariantOperations::remap_data_based_on_genotype(const std::vector<int>& input_data,
    const uint64_t input_call_idx,
    const CombineAllelesLUT& alleles_LUT,
    const unsigned num_merged_alleles, bool NON_REF_exists, const unsigned ploidy,
    RemappedDataWrapperBase& remapped_data,
    std::vector<uint64_t>& num_calls_with_valid_data, int missing_value,
    std::vector<int>& remapped_allele_idx_vec_for_current_gt_combination, std::vector<std::pair<int, int> >& ploidy_index_allele_index_stack,
    std::vector<int>& input_call_allele_idx_vec_for_current_gt_combination);
//Same reason - vcfdiff.cc
template
void VariantOperations::remap_data_based_on_genotype_general(const std::vector<int>& input_data,
    const uint64_t input_call_idx,
    const CombineAllelesLUT& alleles_LUT,
    const unsigned num_merged_alleles, bool NON_REF_exists, const unsigned ploidy,
    RemappedDataWrapperBase& remapped_data,
    std::vector<uint64_t>& num_calls_with_valid_data, int missing_value,
    std::vector<int>& remapped_allele_idx_vec_for_current_gt_combination,
    std::vector<std::pair<int, int> >& ploidy_index_allele_index_stack,
    std::vector<int>& input_call_allele_idx_vec_for_current_gt_combination,
    remap_operator_function_type<int> op
);
