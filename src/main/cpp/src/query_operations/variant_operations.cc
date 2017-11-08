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
#include "query_variants.h"

#ifndef HTSDIR
uint32_t bcf_float_missing    = 0x7F800001;
uint32_t bcf_float_vector_end = 0x7F800002;
#endif
fi_union bcf_float_missing_union = { .i = bcf_float_missing };
fi_union bcf_float_vector_end_union = { .i = bcf_float_vector_end };

//Remapper classes
template<class DataType>
void RemappedMatrix<DataType>::resize(uint64_t num_rows, uint64_t num_columns, DataType init_value)
{
  m_matrix.resize(num_rows);
  for(auto& row : m_matrix)
  {
    row.resize(num_columns);
    for(auto i=0ull;i<num_columns;++i)
      row[i] = init_value;
  }
}

//Each row corresponds to allele/gt idx (all samples), each column corresponds to sample
template<class DataType>
void* RemappedMatrix<DataType>::put_address(uint64_t input_call_idx, unsigned allele_or_gt_idx)
{
  assert(allele_or_gt_idx < m_matrix.size());
  assert(input_call_idx < m_matrix[allele_or_gt_idx].size());
  return reinterpret_cast<void*>(&(m_matrix[allele_or_gt_idx][input_call_idx]));
}

void* RemappedVariant::put_address(uint64_t input_call_idx, unsigned allele_or_gt_idx)
{
  auto& curr_call = m_variant->get_call(input_call_idx);
  assert(curr_call.is_valid());
  auto& field = curr_call.get_field(m_queried_field_idx);
  assert(field.get());  //not null
  return reinterpret_cast<void*>(field->get_address(allele_or_gt_idx)); //returns pointer to the k-th element
}

/*
 * @brief - get the longest reference allele among all variants at this position and store its value in merged_reference_allele
 * For example, if we have the reference alleles T (SNP) and TG (deletion) in two GVCFs at the same location, the reference allele
 * in the merged GVCF should be TG. This modifies the alt alleles (as explained below)
 * Does a sanity check - ref alleles should be a prefix of the merged
 * @param reference_vector - vector of REF strings
 * @param merged_reference_allele - string to store the merged reference
 */
void VariantOperations::merge_reference_allele(const Variant& variant, const VariantQueryConfig& query_config, 
    std::string& merged_reference_allele)
{
  auto* longer_ref = &merged_reference_allele;
  auto merged_ref_length = merged_reference_allele.length();
  if(merged_ref_length == 0u)
  {
    merged_reference_allele = "N";
    merged_ref_length = 1u;
  }
  //assert(variant.get_query_config());
  //const VariantQueryConfig& query_config = *(variant.get_query_config());
  //Iterate over valid calls
  for(const auto& curr_valid_call : variant)
  {
    //If call begins before the current variant begin, its reference is useless
    if(curr_valid_call.get_column_begin() < variant.get_column_begin())
      continue;
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
    if(!CHECK_IN_THE_MIDDLE_REF(merged_reference_allele) && !CHECK_IN_THE_MIDDLE_REF(curr_ref) && longer_ref->find(*shorter_ref) != 0)
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
 * @param merged_alt_alleles vector of merged alt allele strings
 */
void VariantOperations::merge_alt_alleles(const Variant& variant,
    const VariantQueryConfig& query_config,
    const std::string& merged_reference_allele,
    CombineAllelesLUT& alleles_LUT, std::vector<std::string>& merged_alt_alleles, bool& NON_REF_exists) {
  // marking non_reference_allele as already seen will ensure it's not included in the middle
  auto seen_alleles = std::unordered_map<std::string, int>{{g_vcf_NON_REF,-1}};
  merged_alt_alleles.clear();
  auto merged_reference_length = merged_reference_allele.length();
  //invalidate all existing mappings in the LUT
  alleles_LUT.reset_luts();
  //vector to store idx mappings for NON_REF allele, update LUT at end as #ALT alleles are not known till end
  //Set everything to -1 (invalid mapping)
  auto input_non_reference_allele_idx = std::vector<int>(variant.get_num_calls(), -1);
  auto merged_allele_idx = 1u;	//why 1, ref is index 0, alt begins at 1
  NON_REF_exists = false;       //by default, assume NON_REF does not exist
  //Get VariantQueryConfig
  //assert(variant.get_query_config());
  //const VariantQueryConfig& query_config = *(variant.get_query_config());
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
      {
        input_non_reference_allele_idx[curr_call_idx_in_variant] = input_allele_idx;
        NON_REF_exists = true;
      }
      else
      {
        auto* allele_ptr = &allele;
        if(is_suffix_needed && !VariantUtils::is_symbolic_allele(allele))
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
  if(NON_REF_exists)    //if NON_REF allele exists
  {
    // always want non_reference_allele to be last
    merged_alt_alleles.push_back(g_vcf_NON_REF);
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
  }
}

/*
   Remaps GT field
 */
void VariantOperations::remap_GT_field(const std::vector<int>& input_GT, std::vector<int>& output_GT,
    const CombineAllelesLUT& alleles_LUT, const uint64_t input_call_idx, const unsigned num_merged_alleles, const bool NON_REF_exists,
    const FieldLengthDescriptor& length_descriptor)
{
  assert(input_GT.size() == output_GT.size());
  auto should_store_phase_information = length_descriptor.contains_phase_information();
  auto step = should_store_phase_information ? 2u : 1u;
  for(auto i=0u;i<input_GT.size();i+=step)
  {
    if(is_tiledb_missing_value<int>(input_GT[i]) || input_GT[i] == -1 || is_bcf_missing_value<int>(input_GT[i]))
      output_GT[i] = input_GT[i];
    else
    {
      auto output_allele_idx = alleles_LUT.get_merged_idx_for_input(input_call_idx, input_GT[i]);
      if(alleles_LUT.is_missing_value(output_allele_idx))
      {
        if(NON_REF_exists)
        {
          assert(num_merged_alleles >= 2u);
          output_GT[i] = num_merged_alleles-1u;
        }
        else
          output_GT[i] = -1; //missing allele idx
      }
      else
        output_GT[i] = output_allele_idx;
    }
    if(should_store_phase_information && i+1u<input_GT.size())
      output_GT[i+1u] = input_GT[i+1u];
  }
}

// TODO: Implement the genotyping function
void  VariantOperations::do_dummy_genotyping(Variant& variant, std::ostream& output)
{
  assert(variant.get_query_config());
  const VariantQueryConfig& query_config = *(variant.get_query_config());

  for(VariantCall& valid_call : variant)
    modify_reference_if_in_middle(valid_call, query_config, variant.get_column_begin());
  
  std::string merged_reference_allele;
  merged_reference_allele.reserve(10);
  merge_reference_allele(variant, query_config, merged_reference_allele);

  //initialize to number of samples
  CombineAllelesLUT alleles_LUT { static_cast<unsigned>(variant.get_num_calls()) };
  bool NON_REF_exists = false;
  std::vector<std::string> merged_alt_alleles;
  merge_alt_alleles(variant, query_config, merged_reference_allele, alleles_LUT, merged_alt_alleles, NON_REF_exists);

  //Allocate space for remapped PL
  auto num_calls = variant.get_num_calls();
  auto num_merged_alleles = merged_alt_alleles.size() + 1u;     //for REF
  auto num_gts = (num_merged_alleles*(num_merged_alleles+1))/2;
  //Wrapper to store remapped PLs - row corresponds to a single genotype, column to one sample/Call
  RemappedMatrix<int> remapped_PLs;
  remapped_PLs.resize(num_gts, num_calls, bcf_int32_missing);
  std::vector<uint64_t> num_calls_with_valid_data = std::vector<uint64_t>(num_gts, 0ull);

  //Remap PL
  for (auto valid_calls_iter=variant.begin();valid_calls_iter != variant.end();++valid_calls_iter)
  {
    //Not always in sequence, as invalid calls are skipped
    auto curr_call_idx_in_variant = valid_calls_iter.get_call_idx_in_variant();
    auto* PL_field_ptr =
      get_known_field<VariantFieldPrimitiveVectorData<int>, true>(*valid_calls_iter, *(variant.get_query_config()), 
          GVCF_PL_IDX);
    auto* GT_field_ptr =
      get_known_field<VariantFieldPrimitiveVectorData<int>, true>(*valid_calls_iter, *(variant.get_query_config()),
          GVCF_GT_IDX);
    auto length_descriptor = query_config.get_length_descriptor_for_query_attribute_idx(
        query_config.get_query_idx_for_known_field_enum(GVCF_GT_IDX)
        );
    auto ploidy = (GT_field_ptr && GT_field_ptr->is_valid())
      ? length_descriptor.get_ploidy(GT_field_ptr->get().size()) : 2u;
    if(PL_field_ptr && PL_field_ptr->is_valid())
    {
      std::vector<int> tmp_allele_idx_vec;
      std::vector<int> tmp_input_call_allele_idx_vec;
      std::vector<std::pair<int, int> > tmp_stack_vec;
      auto& input_pl_vector = PL_field_ptr->get();
      VariantOperations::remap_data_based_on_genotype<int>(input_pl_vector, curr_call_idx_in_variant,
          alleles_LUT,
          num_merged_alleles, NON_REF_exists, ploidy,
          remapped_PLs,  num_calls_with_valid_data, bcf_int32_missing,
          tmp_allele_idx_vec, tmp_stack_vec,
          tmp_input_call_allele_idx_vec);
    }
  }
  //Compute medians
  std::vector<int> median_vector;
  median_vector.resize(num_gts);
  for(auto i=0u;i<num_gts;++i)
    if(num_calls_with_valid_data[i] == 0ull)
      median_vector[i] = bcf_int32_missing;
    else
    {
      auto& curr_PL_vector = remapped_PLs.get()[i];
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

//VariantOperator functions

//SingleVariantOperatorBase
void SingleVariantOperatorBase::clear()
{
  m_alleles_LUT.reset_luts();
  m_merged_reference_allele.clear();
  for(auto& alt : m_merged_alt_alleles)
    alt.clear();
  m_merged_alt_alleles.clear();
}

void SingleVariantOperatorBase::operate(Variant& variant, const VariantQueryConfig& query_config)
{
  m_merged_reference_allele.resize(0u);
  m_merged_alt_alleles.clear();
  //REF allele
  VariantOperations::merge_reference_allele(variant, query_config, m_merged_reference_allele);
  //ALT alleles
  //set #rows to number of calls
  m_alleles_LUT.resize_luts_if_needed(variant.get_num_calls(), 10u);    //arbitrary non-0 second arg, will be resized correctly anyway
  VariantOperations::merge_alt_alleles(variant, query_config, m_merged_reference_allele, m_alleles_LUT,
      m_merged_alt_alleles, m_NON_REF_exists);
  //is pure reference block if REF is 1 char, and ALT contains only <NON_REF>
  m_is_reference_block_only = (m_merged_reference_allele.length() == 1u && m_merged_alt_alleles.size() == 1u &&
      m_merged_alt_alleles[0] == g_vcf_NON_REF);
  //No remapping is needed if REF is 1 char, and ALT contains only <NON_REF>
  m_remapping_needed = !m_is_reference_block_only;
}

void InterestingLocationsPrinter::operate(Variant& variant, const VariantQueryConfig& query_config)
{
  auto num_valid_calls = 0ull;
  auto num_ref_block_calls = 0ull;
  auto num_begin_at_position = 0ull;
  //Valid calls
  for(const auto& curr_call : variant)
  {
    ++num_valid_calls;
    const auto& ref = get_known_field<VariantFieldString, true>(curr_call, query_config, GVCF_REF_IDX);
    const auto&  alt_field = get_known_field<VariantFieldALTData, true>(curr_call, query_config, GVCF_ALT_IDX);
    if(ref->get().length() == 1u && alt_field->get().size() == 1u
        && alt_field->get()[0u].length() == 1u
        && alt_field->get()[0u][0u] == TILEDB_NON_REF_VARIANT_REPRESENTATION[0u])
      ++num_ref_block_calls;
    if(curr_call.get_column_begin() == variant.get_column_begin())
      ++num_begin_at_position;
  }
  (*m_fptr) << variant.get_column_begin()<<" "<<num_valid_calls<<" "<<num_ref_block_calls<<" "
    <<num_begin_at_position<<"\n";
}

//Dummy genotyping operator
void DummyGenotypingOperator::operate(Variant& variant, const VariantQueryConfig& query_config)
{
  variant.set_query_config(&query_config);
  VariantOperations::do_dummy_genotyping(variant, *m_output_stream);
}

//GA4GHOperator functions
GA4GHOperator::GA4GHOperator(const VariantQueryConfig& query_config, const unsigned max_diploid_alt_alleles_that_can_be_genotyped)
  : SingleVariantOperatorBase()
{
  m_GT_query_idx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
  m_max_diploid_alt_alleles_that_can_be_genotyped = max_diploid_alt_alleles_that_can_be_genotyped;
  m_remapped_fields_query_idxs.clear();
  for(auto query_field_idx=0u;query_field_idx<query_config.get_num_queried_attributes();++query_field_idx)
  {
    //Does the length dependent on number of alleles
    if(query_config.get_length_descriptor_for_query_attribute_idx(query_field_idx).is_length_allele_dependent())
      m_remapped_fields_query_idxs.push_back(query_field_idx);
    //GT field
    if(query_config.get_known_field_enum_for_query_idx(query_field_idx) == GVCF_GT_IDX)
      m_GT_query_idx = query_field_idx;
  }
  m_field_handlers.resize(VARIANT_FIELD_NUM_TYPES);
  m_ploidy.resize(query_config.get_num_rows_in_array());
  for(const auto& ti_enum_pair : g_variant_field_type_index_to_enum)
  {
    unsigned variant_field_type_idx = ti_enum_pair.second;
    assert(variant_field_type_idx < m_field_handlers.size());
    //uninitialized
    assert(m_field_handlers[variant_field_type_idx].get() == 0);
    switch(variant_field_type_idx)
    {
      case VARIANT_FIELD_INT:
        m_field_handlers[variant_field_type_idx] = std::move(std::unique_ptr<VariantFieldHandlerBase>(new VariantFieldHandler<int>())); 
        break;
      case VARIANT_FIELD_UNSIGNED:
        m_field_handlers[variant_field_type_idx] = std::move(std::unique_ptr<VariantFieldHandlerBase>(new VariantFieldHandler<unsigned>())); 
        break;
      case VARIANT_FIELD_INT64_T:
        m_field_handlers[variant_field_type_idx] = std::move(std::unique_ptr<VariantFieldHandlerBase>(new VariantFieldHandler<int64_t>())); 
        break;
      case VARIANT_FIELD_UINT64_T:
        m_field_handlers[variant_field_type_idx] = std::move(std::unique_ptr<VariantFieldHandlerBase>(new VariantFieldHandler<uint64_t>())); 
        break;
      case VARIANT_FIELD_FLOAT:
        m_field_handlers[variant_field_type_idx] = std::move(std::unique_ptr<VariantFieldHandlerBase>(new VariantFieldHandler<float>())); 
        break;                                                                          
      case VARIANT_FIELD_DOUBLE:                                                        
        m_field_handlers[variant_field_type_idx] = std::move(std::unique_ptr<VariantFieldHandlerBase>(new VariantFieldHandler<double>())); 
        break;                                                                          
      case VARIANT_FIELD_STRING:                                                        
        m_field_handlers[variant_field_type_idx] = std::move(std::unique_ptr<VariantFieldHandlerBase>(new VariantFieldHandler<std::string>())); 
        break;                                                                          
      case VARIANT_FIELD_CHAR:                                                          
        m_field_handlers[variant_field_type_idx] = std::move(std::unique_ptr<VariantFieldHandlerBase>(new VariantFieldHandler<char>())); 
        break;
    }
  }
  //Set common fields - REF and ALT for now
  m_remapped_variant.resize_common_fields(2u);
  m_remapped_variant.set_common_field(0u, query_config.get_query_idx_for_known_field_enum(GVCF_REF_IDX), 0);
  m_remapped_variant.set_common_field(1u, query_config.get_query_idx_for_known_field_enum(GVCF_ALT_IDX), 0);
}

std::unique_ptr<VariantFieldHandlerBase>& GA4GHOperator::get_handler_for_type(std::type_index ty)
{
  //Get Enum Idx from VariantFieldTypeEnum
  assert(g_variant_field_type_index_to_enum.find(ty) != g_variant_field_type_index_to_enum.end());
  unsigned variant_field_type_enum = g_variant_field_type_index_to_enum[ty];
  //Check that valid handler exists
  assert(variant_field_type_enum < m_field_handlers.size());
  return m_field_handlers[variant_field_type_enum];
}

void GA4GHOperator::operate(Variant& variant, const VariantQueryConfig& query_config)
{
  //Compute merged REF and ALT
  SingleVariantOperatorBase::operate(variant, query_config);
  //Copy variant to m_remapped_variant - only simple elements, not all fields
  m_remapped_variant.deep_copy_simple_members(variant);
  //Setup code for re-ordering PL/AD etc field elements in m_remapped_variant
  unsigned num_merged_alleles = m_merged_alt_alleles.size()+1u;        //+1 for REF allele
  //Known fields that need to be re-mapped
  if(m_remapping_needed)
  {
    //if GT field is queried
    if(m_GT_query_idx != UNDEFINED_ATTRIBUTE_IDX_VALUE)
    {
      auto GT_length_descriptor = query_config.get_length_descriptor_for_query_attribute_idx(m_GT_query_idx);
      //Valid calls
      for(auto iter=m_remapped_variant.begin();iter!=m_remapped_variant.end();++iter)
      {
        auto& remapped_call = *iter;
        auto curr_call_idx_in_variant = iter.get_call_idx_in_variant();
        m_ploidy[curr_call_idx_in_variant] = 0u;
        auto& remapped_field = remapped_call.get_field(m_GT_query_idx);
        auto& orig_field = variant.get_call(curr_call_idx_in_variant).get_field(m_GT_query_idx);
        copy_field(remapped_field, orig_field);
        if(remapped_field.get() && remapped_field->is_valid())      //Not null
        {
          auto& input_GT =
            variant.get_call(curr_call_idx_in_variant).get_field<VariantFieldPrimitiveVectorData<int>>(m_GT_query_idx)->get();
          auto& output_GT =
            remapped_call.get_field<VariantFieldPrimitiveVectorData<int>>(m_GT_query_idx)->get();
          VariantOperations::remap_GT_field(input_GT, output_GT, m_alleles_LUT, curr_call_idx_in_variant,
              num_merged_alleles, m_NON_REF_exists, GT_length_descriptor);
          m_ploidy[curr_call_idx_in_variant] = GT_length_descriptor.get_ploidy(input_GT.size());
        }
      }
    }
    for(auto query_field_idx : m_remapped_fields_query_idxs)
    {
      auto length_descriptor = query_config.get_length_descriptor_for_query_attribute_idx(query_field_idx);
      //field length depends on #alleles
      assert(length_descriptor.is_length_allele_dependent());
      //Fields such as PL should be skipped, if the #alleles is above a threshold
      if(length_descriptor.is_length_genotype_dependent()
        && too_many_alt_alleles_for_genotype_length_fields(num_merged_alleles-1u))        //#alt = merged-1
      {
        std::cerr << "Column "<<variant.get_column_begin() <<" has too many alleles in the combined VCF record : "<<num_merged_alleles-1
          << " : current limit : "<<m_max_diploid_alt_alleles_that_can_be_genotyped
          << ". Fields, such as  PL, with length equal to the number of genotypes will NOT be added for this location.\n";
        continue;
      }
      //Remapper for m_remapped_variant
      RemappedVariant remapper_variant(m_remapped_variant, query_field_idx); 
      //Iterate over valid calls - m_remapped_variant and variant have same list of valid calls
      for(auto iter=m_remapped_variant.begin();iter!=m_remapped_variant.end();++iter)
      {
        auto& remapped_call = *iter;
        auto curr_call_idx_in_variant = iter.get_call_idx_in_variant();
        auto& remapped_field = remapped_call.get_field(query_field_idx);
        auto& orig_field = variant.get_call(curr_call_idx_in_variant).get_field(query_field_idx);
        copy_field(remapped_field, orig_field);
        if(remapped_field.get() && remapped_field->is_valid())      //Not null
        {
          auto curr_ploidy = m_ploidy[curr_call_idx_in_variant];
          unsigned num_merged_elements =
            length_descriptor.get_num_elements(num_merged_alleles-1u, curr_ploidy, 0u);  //#alt alleles, current ploidy
          remapped_field->resize(num_merged_elements);
          //Get handler for current type
          auto& handler = get_handler_for_type(query_config.get_element_type(query_field_idx));
          assert(handler.get());
          //Call remap function
          handler->remap_vector_data(
              orig_field, curr_call_idx_in_variant,
              m_alleles_LUT, num_merged_alleles, m_NON_REF_exists, curr_ploidy,
              query_config.get_length_descriptor_for_query_attribute_idx(query_field_idx), num_merged_elements, remapper_variant);
        }
      }
    }
  }
  uint64_t offset = 0;
  //Assign REF and ALT common fields
  auto& REF = m_remapped_variant.get_common_field(0u);
  if(REF.get() == 0)
    REF = std::move(std::unique_ptr<VariantFieldString>(new VariantFieldString()));
  auto REF_ptr = dynamic_cast<VariantFieldString*>(REF.get());
  assert(REF_ptr);
  REF_ptr->set_valid(true);
  offset = 0ull;
  REF_ptr->binary_deserialize(&(m_merged_reference_allele[0]), offset, BCF_VL_FIXED, m_merged_reference_allele.length());
  //ALT
  auto& ALT = m_remapped_variant.get_common_field(1u);
  if(ALT.get() == 0)
    ALT = std::move(std::unique_ptr<VariantFieldALTData>(new VariantFieldALTData()));
  auto ALT_ptr =  dynamic_cast<VariantFieldALTData*>(ALT.get()); 
  assert(ALT_ptr);
  ALT_ptr->set_valid(true);
  auto& ALT_vec = ALT_ptr->get();
  ALT_vec.resize(m_merged_alt_alleles.size());
  for(auto i=0u;i<m_merged_alt_alleles.size();++i)
  {
    auto& orig = m_merged_alt_alleles[i];
    auto alt_length = orig.length();
    auto& curr_copy = ALT_vec[i];
    curr_copy.resize(alt_length);
    memcpy(&(curr_copy[0]), &(orig[0]), alt_length*sizeof(char));
  }
}

void GA4GHOperator::copy_back_remapped_fields(Variant& variant) const
{
  if(m_remapping_needed)
  {
    for(auto query_field_idx : m_remapped_fields_query_idxs)
    {
      //Iterate over valid calls - m_remapped_variant and variant have same list of valid calls
      for(auto iter=m_remapped_variant.begin();iter!=m_remapped_variant.end();++iter)
      {
        auto& remapped_call = *iter;
        auto curr_call_idx_in_variant = iter.get_call_idx_in_variant();
        const auto& remapped_field = remapped_call.get_field(query_field_idx);
        auto& orig_field = variant.get_call(curr_call_idx_in_variant).get_field(query_field_idx);
        copy_field(orig_field, remapped_field);
      }
    }
    //if GT field is queried
    if(m_GT_query_idx != UNDEFINED_ATTRIBUTE_IDX_VALUE)
    {
      //Valid calls
      for(auto iter=m_remapped_variant.begin();iter!=m_remapped_variant.end();++iter)
      {
        auto& remapped_call = *iter;
        auto curr_call_idx_in_variant = iter.get_call_idx_in_variant();
        const auto& remapped_field = remapped_call.get_field(m_GT_query_idx);
        auto& orig_field = variant.get_call(curr_call_idx_in_variant).get_field(m_GT_query_idx);
        copy_field(orig_field, remapped_field);
      }
    }
  }
  variant.copy_common_fields(m_remapped_variant);
}


//Single cell operators
ColumnHistogramOperator::ColumnHistogramOperator(uint64_t begin, uint64_t end, uint64_t bin_size)
  : SingleCellOperatorBase()
{
  m_bin_size = bin_size;
  m_begin_column = begin;
  m_end_column = end;
  assert(end >= begin);
  auto num_bins = (end - begin)/bin_size + 1;
  m_bin_counts_vector.resize(num_bins);
  memset(&(m_bin_counts_vector[0]), 0, num_bins*sizeof(uint64_t));
}

void ColumnHistogramOperator::operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema)
{
  auto call_begin = call.get_column_begin();
  auto bin_idx = call_begin <= m_begin_column ? 0ull
    : call_begin >= m_end_column ? m_bin_counts_vector.size()-1 
    : (call_begin - m_begin_column)/m_bin_size;
  assert(bin_idx < m_bin_counts_vector.size());
  ++(m_bin_counts_vector[bin_idx]);
}

void ColumnHistogramOperator::operate_on_columnar_cell(const GenomicsDBColumnarCell& cell, const VariantQueryConfig& query_config,
        const VariantArraySchema& schema)
{
  auto call_begin = static_cast<uint64_t>(cell.get_coordinates()[1]);
  auto bin_idx = call_begin <= m_begin_column ? 0ull
    : call_begin >= m_end_column ? m_bin_counts_vector.size()-1
    : (call_begin - m_begin_column)/m_bin_size;
  assert(bin_idx < m_bin_counts_vector.size());
  ++(m_bin_counts_vector[bin_idx]);
}

bool ColumnHistogramOperator::equi_partition_and_print_bins(uint64_t num_bins, std::ostream& fptr) const
{
  if(num_bins >= m_bin_counts_vector.size())
  {
    std::cerr << "Requested #equi bins is smaller than allocated bin counts vector, returning\n";
    return false;
  }
  auto total_count = 0ull;
  for(auto val : m_bin_counts_vector)
    total_count += val;
  auto count_per_bin = ((double)total_count)/num_bins;
  fptr << "Total "<<total_count<<" #bins "<<num_bins<<" count/bins "<< std::fixed << std::setprecision(1) << count_per_bin <<"\n";
  for(auto i=0ull;i<m_bin_counts_vector.size();)
  {
    auto j = i;
    auto curr_bin_total = 0ull;
    for(;curr_bin_total<count_per_bin && j<m_bin_counts_vector.size();curr_bin_total+=m_bin_counts_vector[j],++j);
    assert(j > i);
    fptr << m_begin_column+i*m_bin_size << "," <<m_begin_column+j*m_bin_size-1 <<"," << curr_bin_total << "\n";
    i = j;
  }
  fptr << "\n";
  return true;
}

void modify_reference_if_in_middle(VariantCall& curr_call, const VariantQueryConfig& query_config, uint64_t current_start_position)
{
  //If the call's column is before the current_start_position, then REF is not valid, set it to "N" (unknown/don't care)
  if(curr_call.get_column_begin() < current_start_position) 
  {
    auto* REF_ptr = get_known_field<VariantFieldString,true>
      (curr_call, query_config, GVCF_REF_IDX);
    REF_ptr->get() = "N";
  }
}

void VariantCallPrintOperator::operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema)
{
  if(m_num_calls_printed > 0ull)
    (*m_fptr) << ",\n";
  call.print(*m_fptr, &query_config, m_indent_prefix, m_vid_mapper);
  ++m_num_calls_printed;
}

void VariantCallPrintOperator::operate_on_columnar_cell(const GenomicsDBColumnarCell& cell, const VariantQueryConfig& query_config,
    const VariantArraySchema& schema)
{
  if(cell.at_new_query_column_interval())
  {
    if(m_num_query_intervals_printed > 0u)
    {
      (*m_fptr) << "\n";
      (*m_fptr) << m_indent_prefix_plus_one << "]\n";
      (*m_fptr) << m_indent_prefix << "},\n";
    }
    m_num_calls_printed = 0u;
    (*m_fptr) << m_indent_prefix << "{\n";
    auto curr_query_column_interval_idx = cell.get_current_query_column_interval_idx();
    auto begin = (query_config.get_num_column_intervals() > 0u)
      ? query_config.get_column_begin(curr_query_column_interval_idx) : 0ll;
    auto end = (query_config.get_num_column_intervals() > 0u)
      ? query_config.get_column_end(curr_query_column_interval_idx) : INT64_MAX-1ll;
    (*m_fptr) << m_indent_prefix_plus_one << "\"query_interval\": [ "<< begin <<", "<< end << " ],\n";
    (*m_fptr) << m_indent_prefix_plus_one << "\"variant_calls\": [\n";
  }
  if(m_num_calls_printed > 0ull)
    (*m_fptr) << ",\n";
  cell.print(*m_fptr, &query_config, m_indent_prefix_plus_two, m_vid_mapper);
  ++m_num_calls_printed;
  ++m_num_query_intervals_printed;
}

void VariantCallPrintOperator::finalize()
{
  if(m_num_query_intervals_printed > 0u)
  {
    (*m_fptr) << "\n";
    (*m_fptr) << m_indent_prefix_plus_one << "]\n";
    (*m_fptr) << m_indent_prefix << "}";
  }
}

void VariantCallPrintCSVOperator::operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema)
{
  auto& fptr = *m_fptr;
  fptr << call.get_row_idx();
  fptr << "," << call.get_column_begin();
  fptr << "," << call.get_column_end();
  //First field is always END - ignore
  for(auto i=1ull;i<query_config.get_num_queried_attributes();++i)
  {
    fptr << ",";
    if(call.get_field(i).get() && call.get_field(i)->is_valid())
    {
      //ALT is handled by concatenating elements with '|' as the separator
      if(query_config.get_known_field_enum_for_query_idx(i) == GVCF_ALT_IDX)
      {
        auto ptr = call.get_field<VariantFieldALTData>(i);
        assert(ptr);
        const auto& alt_vector = ptr->get();
        for(auto j=0u;j<alt_vector.size();++j)
        {
          if(j > 0u)
            fptr << '|';
          fptr << alt_vector[j];
        }
      }
      else
        call.get_field(i)->print_csv(fptr);
    }
    else
    {
      auto schema_idx = query_config.get_schema_idx_for_query_idx(i);
      if(schema.is_variable_length_field(schema_idx))
      {
        auto variant_field_type_enum_idx = VariantFieldTypeUtil::get_variant_field_type_enum_for_variant_field_type(schema.type(schema_idx));
        if(variant_field_type_enum_idx != VariantFieldTypeEnum::VARIANT_FIELD_STRING &&
            variant_field_type_enum_idx != VariantFieldTypeEnum::VARIANT_FIELD_CHAR)
          fptr << "0";
      }
      else
      {
        for(auto j=0;j<schema.val_num(schema_idx)-1;++j)
          fptr << ",";
      }
    }
  }
  fptr << "\n";
}

void VariantCallPrintCSVOperator::operate_on_columnar_cell(const GenomicsDBColumnarCell& cell, const VariantQueryConfig& query_config,
        const VariantArraySchema& schema)
{
  cell.print_csv(*m_fptr, &query_config);
}

//AlleleCountOperator
AlleleCountOperator::AlleleCountOperator(const VidMapper& vid_mapper, const VariantQueryConfig& query_config)
{
  m_vid_mapper = &vid_mapper;
  m_GT_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_GT_IDX);
  if(m_GT_query_idx == UNDEFINED_ATTRIBUTE_IDX_VALUE)
    throw VariantOperationException("GT field must be queried for AlleleCountOperator");
  m_REF_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_REF_IDX);
  if(m_REF_query_idx == UNDEFINED_ATTRIBUTE_IDX_VALUE)
    throw VariantOperationException("REF field must be queried for AlleleCountOperator");
  m_ALT_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_ALT_IDX);
  if(m_ALT_query_idx == UNDEFINED_ATTRIBUTE_IDX_VALUE)
    throw VariantOperationException("ALT field must be queried for AlleleCountOperator");
  auto GT_length_descriptor = query_config.get_length_descriptor_for_query_attribute_idx(m_GT_query_idx);
  assert(GT_length_descriptor.is_length_ploidy_dependent());
  m_GT_step_value = GT_length_descriptor.get_ploidy_step_value();
}

//Unoptimized iterator only traverses single query position/interval
void AlleleCountOperator::operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema)
{
  m_column_to_REF_ALT_to_count_vec.resize(1u);
  auto REF_ptr = call.get_field<VariantFieldString>(m_REF_query_idx);
  auto ALT_ptr = call.get_field<VariantFieldALTData>(m_ALT_query_idx);
  auto GT_ptr = call.get_field<VariantFieldPrimitiveVectorData<int>>(m_GT_query_idx);
  //If any field is invalid, return
  if(!(REF_ptr && REF_ptr->is_valid())
      || !(ALT_ptr && ALT_ptr->is_valid())
      || !(GT_ptr && GT_ptr->is_valid())
      )
    return;
  auto& REF_ALT_to_count_map = get_REF_ALT_to_count_map(call.get_column_begin());
  auto GT_vec = GT_ptr->get();
  auto REF_string = REF_ptr->get();
  auto ALT_vec = ALT_ptr->get();
  //Iterate over GT field
  for(auto i=0u;i<GT_vec.size();i+=m_GT_step_value)
  {
    auto curr_GT_value = GT_vec[i];
    if(is_bcf_valid_value<int>(curr_GT_value) && curr_GT_value > 0) //ignore REF GT
    {
      auto ALT_idx = curr_GT_value-1;
      assert(static_cast<size_t>(ALT_idx) < ALT_vec.size());
      auto REF_ALT_pair = std::pair<std::string, std::string>(
          REF_string, ALT_vec[ALT_idx]);
      normalize_REF_ALT_pair(REF_ALT_pair);
      auto pair_iter = REF_ALT_to_count_map.find(REF_ALT_pair);
      if(pair_iter == REF_ALT_to_count_map.end())
        pair_iter = REF_ALT_to_count_map.insert(std::pair<std::pair<std::string, std::string>, uint64_t>(REF_ALT_pair, 1u)).first;
      else
        ++((*pair_iter).second);
    }
  }
}

void AlleleCountOperator::operate_on_columnar_cell(const GenomicsDBColumnarCell& cell, const VariantQueryConfig& query_config,
        const VariantArraySchema& schema)
{
  if(cell.at_new_query_column_interval())
    m_column_to_REF_ALT_to_count_vec.emplace_back();
  //If any field is invalid, return
  if(!cell.is_valid(m_REF_query_idx) || !cell.is_valid(m_ALT_query_idx)
      || !cell.is_valid(m_GT_query_idx))
    return;
  auto REF_ptr = cell.get_field_ptr_for_query_idx<char>(m_REF_query_idx);
  auto REF_length = cell.get_field_length(m_REF_query_idx);
  auto ALT_ptr = cell.get_field_ptr_for_query_idx<char>(m_ALT_query_idx);
  auto ALT_length = cell.get_field_length(m_ALT_query_idx);
  //ALT is delimited string
  m_cell_ALT_offsets.clear(); //no deallocation
  m_cell_ALT_offsets.push_back(0u); //first ALT allele begins at 0
  auto curr_ALT_ptr = ALT_ptr;
  auto remaining_bytes = ALT_length;
  while(remaining_bytes)
  {
    auto next_ptr = reinterpret_cast<const char*>(memchr(curr_ALT_ptr, TILEDB_ALT_ALLELE_SEPARATOR[0], remaining_bytes));
    auto next_offset = next_ptr
      ? (reinterpret_cast<const char*>(next_ptr)-reinterpret_cast<const char*>(ALT_ptr))+1u //+1 for the delim
      : ALT_length+1u; //+1 for the delim
    //ensures that the last offset element == ALT_length+1
    m_cell_ALT_offsets.push_back(next_offset);
    remaining_bytes = next_ptr ? ALT_length-next_offset : 0;
    curr_ALT_ptr = next_ptr+1u; //skip past the delimiter
  }
  //Find <REF,ALT> -> count map for current column
  auto coords = cell.get_coordinates();
  auto& REF_ALT_to_count_map = get_REF_ALT_to_count_map(coords[1]);
  //Iterate over GT field
  auto GT_ptr = cell.get_field_ptr_for_query_idx<int>(m_GT_query_idx);
  auto GT_length = cell.get_field_length(m_GT_query_idx);
  for(auto i=0u;i<GT_length;i+=m_GT_step_value)
  {
    auto curr_GT_value = GT_ptr[i];
    if(is_bcf_valid_value<int>(curr_GT_value) && curr_GT_value > 0) //ignore REF GT
    {
      auto ALT_idx = curr_GT_value-1;
      assert(static_cast<size_t>(ALT_idx+1) < m_cell_ALT_offsets.size());
      auto REF_ALT_pair = std::move(std::pair<std::string, std::string>(
            std::string(REF_ptr, REF_length),
            std::string(ALT_ptr+m_cell_ALT_offsets[ALT_idx],
              m_cell_ALT_offsets[ALT_idx+1]-m_cell_ALT_offsets[ALT_idx]-1))); //-1 to ignore delim char
      normalize_REF_ALT_pair(REF_ALT_pair);
      auto pair_iter = REF_ALT_to_count_map.find(REF_ALT_pair);
      if(pair_iter == REF_ALT_to_count_map.end())
        pair_iter = REF_ALT_to_count_map.insert(std::pair<std::pair<std::string, std::string>, uint64_t>(REF_ALT_pair, 1u)).first;
      else
        ++((*pair_iter).second);
    }
  }
}

//Normalize ALT before inserting into map
//For example, if REF=TGG ALT=T,AGG  (deletion and SNV), then normalize the SNV to T->A
void AlleleCountOperator::normalize_REF_ALT_pair(std::pair<std::string, std::string>& REF_ALT_pair)
{
  auto REF_length = REF_ALT_pair.first.length();
  auto curr_ALT_length = REF_ALT_pair.second.length();
  auto contains_deletion = (REF_length > 1u);
  if(contains_deletion && curr_ALT_length)
  {
    auto REF_suffix_length = 0u;
    if(VariantUtils::is_symbolic_allele(REF_ALT_pair.second))
      REF_ALT_pair.first.resize(1u); //only store 1 bp for REF in case of symbolic alleles
    else
    {
      if(curr_ALT_length == REF_length) //SNV
      {
        //Last n-1 chars are same
        REF_suffix_length = REF_length-1u;
      }
      else
        if(curr_ALT_length > REF_length) //insertion
        {
          //suffix of ALT and REF must be the same for insertion
          //normalized insertion  A -> ATTG
          //With deletion    ACCG -> ATTGCCG,A
          REF_suffix_length = REF_length-1u;
        }
        else //deletion
        {
          //Could be 2 deletions in the same cell
          //ATCCCG -> ACCCG,A
          //Normalized AT->A and ATCCCG->A
          if(curr_ALT_length > 1u)
            REF_suffix_length = curr_ALT_length-1u;
        }
      assert(REF_ALT_pair.first.substr(REF_length-REF_suffix_length, REF_suffix_length)
          == REF_ALT_pair.second.substr(curr_ALT_length-REF_suffix_length, REF_suffix_length));
      //Chop off suffix
      REF_ALT_pair.first.resize(REF_length-REF_suffix_length);
      REF_ALT_pair.second.resize(curr_ALT_length-REF_suffix_length);
    }
  }
}

std::map<std::pair<std::string, std::string>, uint64_t>& AlleleCountOperator::get_REF_ALT_to_count_map(
    const int64_t curr_column)
{
  auto& column_to_REF_ALT_to_count_map = m_column_to_REF_ALT_to_count_vec.back();
  auto column_to_REF_ALT_to_count_map_iter = column_to_REF_ALT_to_count_map.find(curr_column);
  //If missing, place empty map
  if(column_to_REF_ALT_to_count_map_iter == column_to_REF_ALT_to_count_map.end())
    column_to_REF_ALT_to_count_map_iter = column_to_REF_ALT_to_count_map.insert(std::make_pair(curr_column,
          std::map<std::pair<std::string, std::string>, uint64_t>())).first;
  return (*column_to_REF_ALT_to_count_map_iter).second;
}

void AlleleCountOperator::print_allele_counts(std::ostream& fptr) const
{
  std::string indent_string = "    ";
  std::string curr_indent_string = indent_string;
  //fptr << "{\n";
  //fptr << indent_string << "\"results\": [\n";
  for(const auto& curr_column_to_REF_ALT_to_count : m_column_to_REF_ALT_to_count_vec)
  {
    //curr_indent_string = indent_string + indent_string;
    //fptr << curr_indent_string << "{\n"
    for(const auto& curr_column_to_REF_ALT_to_count_entry : curr_column_to_REF_ALT_to_count)
    {
      for(const auto& curr_REF_ALT_to_count_entry : curr_column_to_REF_ALT_to_count_entry.second)
        fptr << curr_column_to_REF_ALT_to_count_entry.first << " "
          << curr_REF_ALT_to_count_entry.first.first << " "
          << curr_REF_ALT_to_count_entry.first.second << " "
          << curr_REF_ALT_to_count_entry.second <<"\n";
    }
    //fptr << curr_indent_string << "{\n"
  }
  //fptr << indent_string << "]\n"
  //fptr << "}\n";
}
