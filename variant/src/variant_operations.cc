#include "variant_operations.h"

/*
 * @brief - get the longest reference allele among all variants at this position and store its value in merged_reference_allele
 * For example, if we have the reference alleles T (SNP) and TG (deletion) in two GVCFs at the same location, the reference allele
 * in the merged GVCF should be TG. This modifies the alt alleles (as explained below)
 * Does a sanity check - ref alleles should be a prefix of the merged
 * @param reference_vector - vector of REF strings
 * @param merged_reference_allele - string to store the merged reference
 */
void VariantOperations::merge_reference_allele(const Variant& variant, std::string& merged_reference_allele)
{
  auto* longer_ref = &merged_reference_allele;
  auto merged_ref_length = merged_reference_allele.length();
  assert(variant.get_query_config());
  const VariantQueryConfig& query_config = *(variant.get_query_config());
  //Iterate over valid calls
  for(const auto& curr_valid_call : variant)
  {
    auto& curr_ref = get_known_field<VariantFieldString, true>(curr_valid_call, query_config, GVCF_REF_IDX)->get();
    auto curr_ref_length = curr_ref.length();
    auto* shorter_ref = &curr_ref;
    const auto is_curr_ref_longer = (curr_ref_length > merged_ref_length);
    if(is_curr_ref_longer)
    {
      longer_ref = const_cast<std::string*>(&curr_ref);
      shorter_ref = &merged_reference_allele;
    }
#ifdef DEBUG
    //sanity check only - the shorter ref must be a prefix of the longer ref (since they begin at the same location)
    if(longer_ref->find(*shorter_ref) != 0)
    {
      throw std::invalid_argument(std::string{"When combining variants at a given position, the shorter reference allele should be a prefix of the longer reference allele: \'"} + *shorter_ref + " , " + *longer_ref);
    }
#endif
    if(is_curr_ref_longer)
    {
      if(curr_ref_length >= merged_reference_allele.capacity())
	merged_reference_allele.reserve(2*curr_ref_length+1);	//why 2, why not?
      if(merged_ref_length > 0 && CHECK_IN_THE_MIDDLE_REF(merged_reference_allele))
        merged_reference_allele = curr_ref;
      else      //append remaining chars to merged reference
        merged_reference_allele.append(curr_ref, merged_ref_length, curr_ref_length - merged_ref_length);
      merged_ref_length = curr_ref_length;
      longer_ref = &merged_reference_allele;
    }
    else
      if(CHECK_IN_THE_MIDDLE_REF(merged_reference_allele) && !CHECK_IN_THE_MIDDLE_REF(curr_ref))
        merged_reference_allele = curr_ref;
  }
}

/*
 * @brief - remap alt alleles for each input to the longer merged reference
 * For example, if we have the reference alleles T (SNP) and TG (deletion) in two GVCFs at the same location, the reference allele
 * in the merged GVCF should be TG. If the alt alleles were G and T respectively, then the alt alleles in the merged variant
 * become GG,T
 * @param variant - self explanatory
 * @param merged_reference_allele - the merged reference produced by merge_reference_allele
 * @param alleles_LUT -  LUT containing mapping of allele idxs between merged variant and input variants
 * @return vector of merged alt allele strings
 */
const std::vector<std::string>  VariantOperations::merge_alt_alleles(const Variant& variant,
    const std::string& merged_reference_allele,
    CombineAllelesLUT& alleles_LUT) {
  // marking non_reference_allele as already seen will ensure it's not included in the middle
  auto seen_alleles = std::unordered_map<std::string, int>{{g_non_reference_allele,-1}};
  auto merged_alt_alleles = std::vector<std::string>{};
  auto merged_reference_length = merged_reference_allele.length();
  //invalidate all existing mappings in the LUT
  alleles_LUT.reset_luts();
  //vector to store idx mappings for NON_REF allele, update LUT at end as #ALT alleles are not known till end
  //Set everything to -1 (invalid mapping)
  auto input_non_reference_allele_idx = std::vector<int>(variant.get_num_calls(), -1);
  auto merged_allele_idx = 1u;	//why 1, ref is index 0, alt begins at 1
  //Get VariantQueryConfig
  assert(variant.get_query_config());
  const VariantQueryConfig& query_config = *(variant.get_query_config());
  //Iterate over valid calls
  for (auto valid_calls_iter=variant.begin();valid_calls_iter != variant.end();++valid_calls_iter)
  {
    const auto& curr_valid_call = *valid_calls_iter;
    //Not always in sequence, as invalid calls are skipped
    auto curr_call_idx_in_variant = valid_calls_iter.get_call_idx_in_variant();
    const auto& curr_reference =
      get_known_field<VariantFieldString, true>(curr_valid_call, query_config, GVCF_REF_IDX)->get();
    const auto& curr_reference_length = curr_reference.length();
    const auto& curr_allele_vector = 
      get_known_field<VariantFieldALTData, true>(curr_valid_call, query_config, GVCF_ALT_IDX)->get();
    auto is_suffix_needed = false;
    auto suffix_length = 0u;
    if(curr_reference_length < merged_reference_length)
    {
      is_suffix_needed = true;
      suffix_length = merged_reference_length - curr_reference_length;
    }
    //mapping for reference allele 0 -> 0
    alleles_LUT.add_input_merged_idx_pair(curr_call_idx_in_variant, 0, 0);
    auto input_allele_idx = 1u;	//why 1, ref is index 0, alt begins at 1
    //copy of allele if needed
    std::string copy_allele;
    for (const auto& allele : curr_allele_vector)
    {
      if(IS_NON_REF_ALLELE(allele))
        input_non_reference_allele_idx[curr_call_idx_in_variant] = input_allele_idx;
      else
      {
        auto* allele_ptr = &allele;
        if(is_suffix_needed)
        {
          copy_allele = allele;
          copy_allele.append(merged_reference_allele, curr_reference_length, suffix_length);
          allele_ptr = &copy_allele; //allele_ptr now points to a copy of var.alt()[l] (+suffix), hence it's safe to move later
        }
        const auto& iter_pos = seen_alleles.find(*allele_ptr);
        if (iter_pos == seen_alleles.end()) { //allele seen for the first time
          seen_alleles[*allele_ptr] = merged_allele_idx;
          //always check whether LUT is big enough for alleles_LUT (since the #alleles in the merged variant is unknown)
          //Most of the time this function will return quickly (just an if condition check)
          alleles_LUT.resize_luts_if_needed(merged_allele_idx + 1); 
          alleles_LUT.add_input_merged_idx_pair(curr_call_idx_in_variant, input_allele_idx, merged_allele_idx);
          if(is_suffix_needed)
            merged_alt_alleles.push_back(std::move(*allele_ptr)); //allele_ptr points to a copy - use move
          else
            merged_alt_alleles.push_back(*allele_ptr);	//allele_ptr points to curr_allele_vector[input_allele_idx] - copy to vector
          ++merged_allele_idx;
        }
        else
          alleles_LUT.add_input_merged_idx_pair(curr_call_idx_in_variant, input_allele_idx, (*iter_pos).second);
      }
      ++input_allele_idx;
    }
  }
  // always want non_reference_allele to be last
  merged_alt_alleles.push_back(g_non_reference_allele);
  auto non_reference_allele_idx = merged_alt_alleles.size(); //why not -1, include reference allele also
  //always check whether LUT is big enough for alleles_LUT (since the #alleles in the merged variant is unknown)
  alleles_LUT.resize_luts_if_needed(non_reference_allele_idx + 1); 
  //Add mappings for non_ref allele
  //Iterate over valid calls
  for (auto valid_calls_iter=variant.begin();valid_calls_iter != variant.end();++valid_calls_iter)
  {
    //Not always in sequence, as invalid calls are skipped
    auto curr_call_idx_in_variant = valid_calls_iter.get_call_idx_in_variant();
    if(input_non_reference_allele_idx[curr_call_idx_in_variant] >= 0)
      alleles_LUT.add_input_merged_idx_pair(curr_call_idx_in_variant, input_non_reference_allele_idx[curr_call_idx_in_variant],
          non_reference_allele_idx);
  }
  return merged_alt_alleles;
}

/*
  Copied from defunct gamgee library
  @param pls - vector of PL values for a given sample as stored in TileDB
  @param input_sample_idx 
  @param alleles_LUT LUT mapping alleles from input to merged alleles list
  @param num_merged_alleles  
  @param remapped_pls - matrix of PLs, each row corresponds to a single genotype, each column corresponds to a sample
  @num_calls_with_valid_data - keeps track of how many samples had valid values for given genotype idx
 */
void  VariantOperations::remap_pl(const std::vector<int32_t>& pls, const uint32_t input_sample_idx,
    const CombineAllelesLUT& alleles_LUT, const unsigned num_merged_alleles,
    std::vector<std::vector<int>>& remapped_pls,
    std::vector<uint64_t>& num_calls_with_valid_data) {
  //index of NON_REF in merged variant
  const auto merged_non_reference_allele_idx = static_cast<int>(num_merged_alleles-1);
  //index of NON_REF in input sample
  const auto input_non_reference_allele_idx = alleles_LUT.get_input_idx_for_merged(input_sample_idx, merged_non_reference_allele_idx);

  for (auto allele_j = 0u; allele_j < num_merged_alleles; ++allele_j) {
    auto input_j_allele = alleles_LUT.get_input_idx_for_merged(input_sample_idx, allele_j);
    if (CombineAllelesLUT::is_missing_value(input_j_allele))	//no mapping found for current allele in input gvcf
    {
      if(CombineAllelesLUT::is_missing_value(input_non_reference_allele_idx))	//input did not have NON_REF allele
      {
	//fill in missing values for all genotypes with allele_j as one component
	for(auto allele_k = allele_j; allele_k < num_merged_alleles;++allele_k)
        {
          auto gt_idx = bcf_alleles2gt(allele_j, allele_k);
	  remapped_pls[gt_idx][input_sample_idx] = bcf_int32_missing;
        }
	continue;	//skip to next value of allele_j
      }
      else //input contains NON_REF allele, use its idx
	input_j_allele = input_non_reference_allele_idx;
    }
    for (auto allele_k = allele_j; allele_k < num_merged_alleles; ++allele_k) {
      auto gt_idx = bcf_alleles2gt(allele_j, allele_k);
      auto input_k_allele = alleles_LUT.get_input_idx_for_merged(input_sample_idx, allele_k);
      if (CombineAllelesLUT::is_missing_value(input_k_allele))	//no mapping found for current allele in input gvcf
      {
	if(CombineAllelesLUT::is_missing_value(input_non_reference_allele_idx))	//input did not have NON_REF allele
	{
	  remapped_pls[gt_idx][input_sample_idx] = bcf_int32_missing; //put missing value
	  continue;	//skip to next value of allele_k
	}
	else //input has NON_REF, use its idx
	  input_k_allele = input_non_reference_allele_idx;
      }
      remapped_pls[gt_idx][input_sample_idx] = pls[bcf_alleles2gt(input_j_allele, input_k_allele)];
      ++(num_calls_with_valid_data[gt_idx]);
    }
  }
}

// TODO: Implement the genotyping function
void  VariantOperations::do_dummy_genotyping(Variant& variant, std::ostream& output)
{
  assert(variant.get_query_config());
  std::string merged_reference_allele;
  merged_reference_allele.reserve(10);
  merge_reference_allele(variant, merged_reference_allele);

  //initialize to number of samples
  CombineAllelesLUT alleles_LUT { static_cast<unsigned>(variant.get_num_calls()) };
  auto& merged_alt_alleles = merge_alt_alleles(variant, merged_reference_allele, alleles_LUT);

  //Allocate space for remapped PL
  auto num_calls = variant.get_num_calls();
  auto num_merged_alleles = merged_alt_alleles.size() + 1u;
  auto num_gts = (num_merged_alleles*(num_merged_alleles+1))/2;
  std::vector<std::vector<int>> remapped_PLs = std::vector<std::vector<int>>(num_gts,
      std::vector<int>(num_calls, bcf_int32_missing));
  std::vector<uint64_t> num_calls_with_valid_data = std::vector<uint64_t>(num_gts, 0ull);

  //Remap PL
  for (auto valid_calls_iter=variant.begin();valid_calls_iter != variant.end();++valid_calls_iter)
  {
    //Not always in sequence, as invalid calls are skipped
    auto curr_call_idx_in_variant = valid_calls_iter.get_call_idx_in_variant();
    auto& input_pl_vector = 
      get_known_field<VariantFieldPrimitiveVectorData<int>, true>(*valid_calls_iter, *(variant.get_query_config()), 
          GVCF_PL_IDX)->get();
    remap_pl(input_pl_vector, curr_call_idx_in_variant,
        alleles_LUT, num_merged_alleles,
        remapped_PLs,  num_calls_with_valid_data);
  }
  //Compute medians
  std::vector<int> median_vector;
  median_vector.resize(num_gts);
  for(auto i=0u;i<num_gts;++i)
    if(num_calls_with_valid_data[i] == 0ull)
      median_vector[i] = bcf_int32_missing;
    else
    {
      auto& curr_PL_vector = remapped_PLs[i];
      auto dec_order_median_idx = (num_calls_with_valid_data[i])/2;
      //auto inc_order_median_idx = num_calls_with_valid_data[i])/2;
      std::nth_element(curr_PL_vector.begin(), curr_PL_vector.begin() + dec_order_median_idx, curr_PL_vector.end(), std::greater<int>());
      //std::nth_element(curr_PL_vector.begin(), curr_PL_vector.begin() + inc_order_median_idx, curr_PL_vector.end());
      median_vector[i] = curr_PL_vector[dec_order_median_idx];
      //median_vector[i] = curr_PL_vector[inc_order_median_idx];
      assert(median_vector[i] != bcf_int32_missing);
    }
  output << variant.get_column_begin() << ",";
  output << merged_reference_allele;
  for(const auto& curr_alt_allele : merged_alt_alleles)
    output << "," << curr_alt_allele;
  for(auto value : median_vector)
    output << "," << value;
  output << "\n";
  return;
}
