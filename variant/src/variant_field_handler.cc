#include "variant_operations.h"

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
      if(CombineAllelesLUT::is_missing_value(input_non_reference_allele_idx))	//input did not have NON_REF allele
      {
        *(reinterpret_cast<DataType*>(remapped_data.put_address(input_call_idx, j))) = (missing_value);
        continue;
      }
      else //input contains NON_REF allele, use its idx
        input_j_allele = input_non_reference_allele_idx;
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
  if(KnownFieldInfo::is_length_genotype_dependent(length_descriptor)) 
    VariantOperations::remap_data_based_on_genotype<DataType>( 
        orig_vector_field_ptr->get(), curr_call_idx_in_variant, 
        alleles_LUT, num_merged_alleles, non_ref_exists, 
        remapper_variant, m_num_calls_with_valid_data, m_bcf_missing_value); 
  else 
    VariantOperations::remap_data_based_on_alleles<DataType>( 
        orig_vector_field_ptr->get(), curr_call_idx_in_variant, 
        alleles_LUT, num_merged_alleles, non_ref_exists, 
        KnownFieldInfo::is_length_only_ALT_alleles_dependent(length_descriptor), 
        remapper_variant, m_num_calls_with_valid_data, m_bcf_missing_value);
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
