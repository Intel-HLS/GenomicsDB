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
  //assert(variant.get_query_config());
  //const VariantQueryConfig& query_config = *(variant.get_query_config());
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
 * @return vector of merged alt allele strings
 */
const std::vector<std::string>  VariantOperations::merge_alt_alleles(const Variant& variant,
    const VariantQueryConfig& query_config,
    const std::string& merged_reference_allele,
    CombineAllelesLUT& alleles_LUT, bool& NON_REF_exists) {
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
  if(NON_REF_exists)    //if NON_REF allele exists
  {
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
  }
  return merged_alt_alleles;
}

/*
   Remaps GT field
 */
void VariantOperations::remap_GT_field(const std::vector<int>& input_GT, std::vector<int>& output_GT,
    const CombineAllelesLUT& alleles_LUT, const uint64_t input_call_idx)
{
  assert(input_GT.size() == output_GT.size());
  for(auto i=0;i<input_GT.size();++i)
  {
    auto output_allele_idx = alleles_LUT.get_merged_idx_for_input(input_call_idx, input_GT[i]);
    assert(!alleles_LUT.is_missing_value(output_allele_idx));
    output_GT[i] = output_allele_idx;
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
  auto& merged_alt_alleles = merge_alt_alleles(variant, query_config, merged_reference_allele, alleles_LUT, NON_REF_exists);

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
    if(PL_field_ptr && PL_field_ptr->is_valid())
    {
      auto& input_pl_vector = PL_field_ptr->get();
      remap_data_based_on_genotype<int>(input_pl_vector, curr_call_idx_in_variant,
          alleles_LUT, num_merged_alleles, NON_REF_exists,
          remapped_PLs,  num_calls_with_valid_data, bcf_int32_missing);
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
  //REF allele
  VariantOperations::merge_reference_allele(variant, query_config, m_merged_reference_allele);
  //ALT alleles
  //set #rows to number of calls
  m_alleles_LUT.resize_luts_if_needed(variant.get_num_calls(), 10u);    //arbitrary non-0 second arg, will be resized correctly anyway
  auto& merged_alt_alleles =  VariantOperations::merge_alt_alleles(variant, query_config, m_merged_reference_allele, m_alleles_LUT,
      m_NON_REF_exists);
  //Move to class member
  m_merged_alt_alleles = std::move(merged_alt_alleles);
}

//Dummy genotyping operator
void DummyGenotypingOperator::operate(Variant& variant, const VariantQueryConfig& query_config)
{
  variant.set_query_config(&query_config);
  VariantOperations::do_dummy_genotyping(variant, *m_output_stream);
}

//GA4GHOperator functions
GA4GHOperator::GA4GHOperator() : SingleVariantOperatorBase()
{
  m_field_handlers.resize(VARIANT_FIELD_NUM_TYPES);
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
}

inline std::unique_ptr<VariantFieldHandlerBase>& GA4GHOperator::get_handler_for_type(std::type_index ty)
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
  //Store copy of variant in vector<Variant>
  m_variants.emplace_back();
  Variant& copy = m_variants[m_variants.size()-1u];
  copy.copy_from_variant(variant);      //Create copy to store altered PL fields
  //Setup code for re-ordering PL/AD etc field elements in copy
  unsigned num_merged_alleles = m_merged_alt_alleles.size()+1u;        //+1 for REF allele
  unsigned num_genotypes = (num_merged_alleles*(num_merged_alleles+1u))/2u;
  for(auto query_field_idx=0u;query_field_idx<query_config.get_num_queried_attributes();++query_field_idx)
  {
    //is known field?
    if(query_config.is_defined_known_field_enum_for_query_idx(query_field_idx))
    {
      const auto* info_ptr = query_config.get_info_for_query_idx(query_field_idx);
      //known field whose length is dependent on #alleles
      if(info_ptr && info_ptr->is_length_allele_dependent())
      {
        unsigned num_merged_elements = info_ptr->get_num_elements_for_known_field_enum(num_merged_alleles-1u, 0u);     //#alt alleles
        //Remapper for copy
        RemappedVariant remapper_variant(copy, query_field_idx); 
        //Iterate over valid calls - copy and variant have same list of valid calls
        for(auto iter=copy.begin();iter!=copy.end();++iter)
        {
          auto& curr_call = *iter;
          auto curr_call_idx_in_variant = iter.get_call_idx_in_variant();
          auto& curr_field = curr_call.get_field(query_field_idx);
          if(curr_field.get() && curr_field->is_valid())      //Not null
          {
            curr_field->resize(num_merged_elements);
            //Get handler for current type
            auto& handler = get_handler_for_type(curr_field->get_element_type());
            assert(handler.get());
            //Call remap function
            handler->remap_vector_data(
                variant.get_call(curr_call_idx_in_variant).get_field(query_field_idx), curr_call_idx_in_variant,
                m_alleles_LUT, num_merged_alleles, m_NON_REF_exists,
                info_ptr->get_length_descriptor(), num_merged_elements, remapper_variant);
          }
        }
      }
      //GT field
      if(query_config.get_known_field_enum_for_query_idx(query_field_idx) == GVCF_GT_IDX)
      {
        //Valid calls
        for(auto iter=copy.begin();iter!=copy.end();++iter)
        {
          auto& curr_call = *iter;
          auto curr_call_idx_in_variant = iter.get_call_idx_in_variant();
          auto& curr_field = curr_call.get_field(query_field_idx);
          if(curr_field.get() && curr_field->is_valid())      //Not null
          {
            auto& input_GT =
              variant.get_call(curr_call_idx_in_variant).get_field<VariantFieldPrimitiveVectorData<int>>(query_field_idx)->get();
            auto& output_GT = 
              curr_call.get_field<VariantFieldPrimitiveVectorData<int>>(query_field_idx)->get();
            VariantOperations::remap_GT_field(input_GT, output_GT, m_alleles_LUT, curr_call_idx_in_variant);
          }
        }
      }
    }
  }
  //Set common fields - REF and ALT for now
  copy.resize_common_fields(2u);
  auto* REF_ptr = new VariantFieldString();
  REF_ptr->set_valid(true);
  REF_ptr->get() = std::move(m_merged_reference_allele);        //get returns string&
  copy.set_common_field(0u, query_config.get_query_idx_for_known_field_enum(GVCF_REF_IDX), REF_ptr);
  //ALT
  auto* ALT_ptr = new VariantFieldALTData();
  ALT_ptr->set_valid(true);
  ALT_ptr->get() = std::move(m_merged_alt_alleles);        //get returns vector<string>&
  copy.set_common_field(1u, query_config.get_query_idx_for_known_field_enum(GVCF_ALT_IDX), ALT_ptr);
  //Do not use m_merged_alt_alleles and m_merged_reference_allele after this point
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

