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

#ifdef HTSDIR

#include "broad_combined_gvcf.h"

//INFO fields
#define MAKE_BCF_INFO_TUPLE(enum_idx, query_idx, variant_type_enum, bcf_type, vcf_field_name, VCF_field_combine_operation) \
  INFO_tuple_type(enum_idx, query_idx, variant_type_enum, bcf_type, vcf_field_name, VCF_field_combine_operation)
#define BCF_INFO_GET_KNOWN_FIELD_ENUM(X) (std::get<0>(X))
#define BCF_INFO_GET_QUERY_FIELD_IDX(X) (std::get<1>(X))
#define BCF_INFO_GET_VARIANT_FIELD_TYPE_ENUM(X) (std::get<2>(X))
#define BCF_INFO_GET_BCF_HT_TYPE(X) (std::get<3>(X))
#define BCF_INFO_GET_VCF_FIELD_NAME(X) (std::get<4>(X))
#define BCF_INFO_GET_VCF_FIELD_COMBINE_OPERATION(X) (std::get<5>(X))
//FORMAT fields
#define MAKE_BCF_FORMAT_TUPLE(enum_idx, query_idx, variant_type_enum, bcf_type, vcf_field_name) \
  FORMAT_tuple_type(enum_idx, query_idx, variant_type_enum, bcf_type, vcf_field_name)
#define BCF_FORMAT_GET_KNOWN_FIELD_ENUM(X) (std::get<0>(X))
#define BCF_FORMAT_GET_QUERY_FIELD_IDX(X) (std::get<1>(X))
#define BCF_FORMAT_GET_VARIANT_FIELD_TYPE_ENUM(X) (std::get<2>(X))
#define BCF_FORMAT_GET_BCF_HT_TYPE(X) (std::get<3>(X))
#define BCF_FORMAT_GET_VCF_FIELD_NAME(X) (std::get<4>(X))

//Static member
const std::unordered_set<char> BroadCombinedGVCFOperator::m_legal_bases({'A', 'T', 'G', 'C'});

//Utility functions for encoding GT field
template<bool phase_information_in_TileDB, bool produce_GT_field>
int encode_GT_element(const int value, const bool is_phased);

//Specialization - phase information exists in TileDB and produce_GT_field
template<>
inline int encode_GT_element<true, true>(const int value, const bool is_phased)
{
  return is_bcf_valid_value(value)
    ? (is_phased ? bcf_gt_phased(value) : bcf_gt_unphased(value))
    : value;
}

//Specialization - phase information exists in TileDB and !produce_GT_field
template<>
inline int encode_GT_element<true, false>(const int value, const bool is_phased)
{
  return is_bcf_valid_value(value)
    ? (is_phased ? (bcf_gt_missing | 1) : bcf_gt_missing)
    : value;
}
//Specialization - phase information doesn't exist in TileDB and produce_GT_field
template<>
inline int encode_GT_element<false, true>(const int value, const bool is_phased)
{
  return is_bcf_valid_value(value) ? bcf_gt_unphased(value) : value;
}

//Specialization - phase information doesn't exist in TileDB and !produce_GT_field
template<>
inline int encode_GT_element<false, false>(const int value, const bool is_phased)
{
  return is_bcf_valid_value(value) ? bcf_gt_missing : value;
}

template<bool phase_information_in_TileDB, bool produce_GT_field>
void encode_GT_vector(int* inout_vec, const uint64_t input_offset,
    const unsigned num_elements_per_sample, uint64_t& output_idx);

//Specialization - phase information exists in TileDB and produce_GT_field
template<>
inline void encode_GT_vector<true, true>(int* inout_vec, const uint64_t input_offset,
    const unsigned num_elements_per_sample, uint64_t& output_idx)
{
  if(num_elements_per_sample > 0u)
    inout_vec[output_idx++] = encode_GT_element<false, true>(inout_vec[input_offset], false);
  //skip over phasing elements
  for(auto k=2u;k<num_elements_per_sample;k+=2u)
    inout_vec[output_idx++] = encode_GT_element<true, true>(inout_vec[input_offset+k], (inout_vec[input_offset+k-1u] > 0));
}

//Specialization - phase information exists in TileDB and !produce_GT_field
template<>
inline void encode_GT_vector<true, false>(int* inout_vec, const uint64_t input_offset,
    const unsigned num_elements_per_sample, uint64_t& output_idx)
{
  if(num_elements_per_sample > 0u)
    inout_vec[output_idx++] = encode_GT_element<false, false>(inout_vec[input_offset], false);
  //skip over phasing elements
  for(auto k=2u;k<num_elements_per_sample;k+=2u)
    inout_vec[output_idx++] = encode_GT_element<true, false>(inout_vec[input_offset+k], (inout_vec[input_offset+k-1u] > 0));
}

//Specialization - phase information doesn't exist in TileDB and produce_GT_field
template<>
inline void encode_GT_vector<false, true>(int* inout_vec, const uint64_t input_offset,
    const unsigned num_elements_per_sample, uint64_t& output_idx)
{
  if(num_elements_per_sample > 0u)
    inout_vec[output_idx++] = encode_GT_element<false, true>(inout_vec[input_offset], false);
  //skip over phasing elements
  for(auto k=1u;k<num_elements_per_sample;++k)
    inout_vec[output_idx++] = encode_GT_element<false, true>(inout_vec[input_offset+k], false);
}

//Specialization - phase information doesn't exist in TileDB and !produce_GT_field
template<>
inline void encode_GT_vector<false, false>(int* inout_vec, const uint64_t input_offset,
    const unsigned num_elements_per_sample, uint64_t& output_idx)
{
  if(num_elements_per_sample > 0u)
    inout_vec[output_idx++] = encode_GT_element<false, false>(inout_vec[input_offset], false);
  //skip over phasing elements
  for(auto k=1u;k<num_elements_per_sample;++k)
    inout_vec[output_idx++] = encode_GT_element<false, false>(inout_vec[input_offset+k], false);
}

BroadCombinedGVCFOperator::BroadCombinedGVCFOperator(VCFAdapter& vcf_adapter, const VidMapper& id_mapper,
    const VariantQueryConfig& query_config,
    const unsigned max_diploid_alt_alleles_that_can_be_genotyped, const bool use_missing_values_only_not_vector_end)
: GA4GHOperator(query_config, max_diploid_alt_alleles_that_can_be_genotyped)
{
  clear();
  if(!id_mapper.is_initialized())
    throw BroadCombinedGVCFException("Id mapper is not initialized");
  if(!id_mapper.is_callset_mapping_initialized())
    throw BroadCombinedGVCFException("Callset mapping in id mapper is not initialized");
  m_query_config = &query_config;
  //Initialize VCF structs
  m_vcf_adapter = &vcf_adapter;
  m_vid_mapper = &id_mapper;
  m_use_missing_values_not_vector_end = use_missing_values_only_not_vector_end;
  m_vcf_hdr = vcf_adapter.get_vcf_header();
  m_bcf_out = bcf_init();
  //vector of char*, to avoid frequent reallocs()
  m_alleles_pointer_buffer.resize(100u);
  //DP INFO field - handle after all FORMAT fields have been processed
  FORMAT_tuple_type DP_INFO_as_FORMAT_tuple;
  auto is_DP_INFO_queried = false;
  //Determine queried INFO and FORMAT fields
  for(auto i=0u;i<query_config.get_num_queried_attributes();++i)
  {
    auto* field_info = m_vid_mapper->get_field_info(query_config.get_query_attribute_name(i));
    if(field_info)
    {
      auto known_field_enum = query_config.is_defined_known_field_enum_for_query_idx(i) ? query_config.get_known_field_enum_for_query_idx(i)
        : UNDEFINED_ATTRIBUTE_IDX_VALUE;
      auto VCF_field_combine_operation = query_config.get_VCF_field_combine_operation_for_query_attribute_idx(i);
      auto add_to_INFO_vector = (field_info->m_is_vcf_INFO_field && known_field_enum != GVCF_END_IDX
          && (known_field_enum != GVCF_DP_IDX || VCF_field_combine_operation != VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_DP) //not DP or not combined as GATK Combine GVCF DP
          && VCF_field_combine_operation != VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MOVE_TO_FORMAT  //not moved to FORMAT
          );
      auto add_to_FORMAT_vector = (field_info->m_is_vcf_FORMAT_field ||
          (field_info->m_is_vcf_INFO_field
           && ((known_field_enum == GVCF_DP_IDX && VCF_field_combine_operation == VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_DP)
             || (VCF_field_combine_operation == VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MOVE_TO_FORMAT)
           )
           )
          );
      if(add_to_INFO_vector)
      {
        if(VCF_field_combine_operation == VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_UNKNOWN_OPERATION)
          std::cerr << "WARNING: No valid combination operation found for INFO field "<<field_info->m_vcf_name<<" - the field will NOT be part of INFO fields in the generated VCF records\n";
        else
        {
          m_INFO_fields_vec.emplace_back(MAKE_BCF_INFO_TUPLE(known_field_enum, i,
                VariantFieldTypeUtil::get_variant_field_type_enum_for_variant_field_type(field_info->m_type_index),
                VariantFieldTypeUtil::get_vcf_field_type_enum_for_variant_field_type(field_info->m_type_index),
                field_info->m_vcf_name,
                VCF_field_combine_operation));
          VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, &id_mapper, field_info->m_vcf_name, BCF_HL_INFO);
        }
      }
      if(add_to_FORMAT_vector)
      {
        auto format_tuple = MAKE_BCF_FORMAT_TUPLE(known_field_enum, i,
            VariantFieldTypeUtil::get_variant_field_type_enum_for_variant_field_type(field_info->m_type_index),
            VariantFieldTypeUtil::get_vcf_field_type_enum_for_variant_field_type(field_info->m_type_index),
            field_info->m_vcf_name);
        if((field_info->m_is_vcf_FORMAT_field)
            || (VCF_field_combine_operation == VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MOVE_TO_FORMAT)
          )
        {
          m_FORMAT_fields_vec.emplace_back(format_tuple);
          VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, &id_mapper, field_info->m_vcf_name, BCF_HL_FMT);
        }
        else //DP INFO
        {
          DP_INFO_as_FORMAT_tuple = std::move(format_tuple);
          is_DP_INFO_queried = true;
          VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, &id_mapper, "DP", BCF_HL_INFO);
        }
      }
    }
  }
  //is FILTER field queried?
  if(query_config.is_defined_query_idx_for_known_field_enum(GVCF_FILTER_IDX))
  {
    for(auto i=0u;i<m_vid_mapper->get_num_fields();++i)
    {
      auto& field_info = m_vid_mapper->get_field_info(i);
      if(field_info.m_is_vcf_FILTER_field)
        VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, m_vid_mapper, field_info.m_vcf_name, BCF_HL_FLT);
    }
  }
  //If DP is queried, it should always be the last field in m_FORMAT_fields_vec
  if(is_DP_INFO_queried)
  {
    //Move to the last element
    m_FORMAT_fields_vec.emplace_back(DP_INFO_as_FORMAT_tuple);
    if(BCF_FORMAT_GET_KNOWN_FIELD_ENUM(m_FORMAT_fields_vec[m_FORMAT_fields_vec.size()-1u]) != GVCF_DP_IDX)
      throw BroadCombinedGVCFException("Last queried FORMAT field should be DP, instead it is "
          +g_known_variant_field_names[BCF_FORMAT_GET_KNOWN_FIELD_ENUM(m_FORMAT_fields_vec[m_FORMAT_fields_vec.size()-1u])]);
  }
  //Qual combine operation
  m_vcf_qual_tuple = MAKE_BCF_INFO_TUPLE(KnownVariantFieldsEnum::GVCF_QUAL_IDX, UNDEFINED_ATTRIBUTE_IDX_VALUE,
                VariantFieldTypeEnum::VARIANT_FIELD_FLOAT,
                BCF_HT_REAL,
                g_known_variant_field_names[KnownVariantFieldsEnum::GVCF_QUAL_IDX],
                VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_UNKNOWN_OPERATION);
  if(query_config.is_defined_query_idx_for_known_field_enum(KnownVariantFieldsEnum::GVCF_QUAL_IDX))
  {
    auto query_field_idx = query_config.get_query_idx_for_known_field_enum(KnownVariantFieldsEnum::GVCF_QUAL_IDX);
    std::get<1>(m_vcf_qual_tuple) = query_field_idx;
    std::get<5>(m_vcf_qual_tuple) = query_config.get_VCF_field_combine_operation_for_query_attribute_idx(query_field_idx);
  }
  //Add missing contig names to template header
  for(auto i=0u;i<m_vid_mapper->get_num_contigs();++i)
  {
    auto& contig_info = m_vid_mapper->get_contig_info(i);
    auto& contig_name = contig_info.m_name;
    if(bcf_hdr_name2id(m_vcf_hdr, contig_name.c_str()) < 0)
    {
      std::string contig_vcf_line = std::string("##contig=<ID=")+contig_name+",length="
        +std::to_string(contig_info.m_length)+">";
      int line_length = 0;
      auto hrec = bcf_hdr_parse_line(m_vcf_hdr, contig_vcf_line.c_str(), &line_length);
      bcf_hdr_add_hrec(m_vcf_hdr, hrec);
      bcf_hdr_sync(m_vcf_hdr);
    }
  }
  //Get contig info for position 0, store curr contig in next_contig and call switch_contig function to do all the setup
  auto curr_contig_flag = m_vid_mapper->get_next_contig_location(-1ll, m_next_contig_name, m_next_contig_begin_position);
  assert(curr_contig_flag);
  switch_contig();
  //Add samples to template header
  std::string callset_name;
  for(auto i=0ull;i<query_config.get_num_rows_to_query();++i)
  {
    auto row_idx = query_config.get_array_row_idx_for_query_row_idx(i);
    auto status = m_vid_mapper->get_callset_name(row_idx, callset_name);
    if(!status || callset_name.empty())
      throw BroadCombinedGVCFException(std::string("No sample/CallSet name specified in JSON file/Protobuf object for TileDB row ")
          + std::to_string(row_idx));
    auto add_sample_status = bcf_hdr_add_sample(m_vcf_hdr, callset_name.c_str());
    if(add_sample_status < 0)
      throw BroadCombinedGVCFException(std::string("Could not add sample ")
          +callset_name+" to the combined VCF/gVCF header");
  }
  bcf_hdr_sync(m_vcf_hdr);
  //Map from vid mapper field idx to hdr field idx
  m_global_field_idx_to_hdr_idx.resize(m_vid_mapper->get_num_fields(), -1);
  for(auto i=0u;i<m_vid_mapper->get_num_fields();++i)
  {
    auto& field_info = m_vid_mapper->get_field_info(i);
    //Could be -1
    m_global_field_idx_to_hdr_idx[i] = bcf_hdr_id2int(m_vcf_hdr, BCF_DT_ID, field_info.m_vcf_name.c_str());
  }
  m_vcf_adapter->print_header();
  //vector of field pointers used for handling remapped fields when dealing with spanning deletions
  //Individual pointers will be allocated later
  m_spanning_deletions_remapped_fields.resize(m_remapped_fields_query_idxs.size());
  //Initialize GT encoding function pointer
  const auto& GT_length_descriptor = query_config.get_length_descriptor_for_query_attribute_idx(m_GT_query_idx);
  if(GT_length_descriptor.contains_phase_information())
  {
    //GATK CombineGVCF does not produce GT field by default - option to produce GT
    if(m_vcf_adapter->produce_GT_field())
      m_encode_GT_vector_function_ptr = encode_GT_vector<true, true>;
    else
      m_encode_GT_vector_function_ptr = encode_GT_vector<true, false>;
  }
  else
  {
    //GATK CombineGVCF does not produce GT field by default - option to produce GT
    if(m_vcf_adapter->produce_GT_field())
      m_encode_GT_vector_function_ptr = encode_GT_vector<false, true>;
    else
      m_encode_GT_vector_function_ptr = encode_GT_vector<false, false>;
  }
}

void BroadCombinedGVCFOperator::clear()
{
  m_curr_contig_name.clear();
  m_next_contig_name.clear();
  m_alleles_pointer_buffer.clear();
  m_INFO_fields_vec.clear();
  m_FORMAT_fields_vec.clear();
  m_MIN_DP_vector.clear();
  m_DP_FORMAT_vector.clear();
  m_spanning_deletions_remapped_fields.clear();
  m_spanning_deletion_remapped_GT.clear();
  m_spanning_deletion_current_genotype.clear();
  m_global_field_idx_to_hdr_idx.clear();
  m_FILTER_idx_vec.clear();
}

bool BroadCombinedGVCFOperator::handle_VCF_field_combine_operation(const Variant& variant,
    const INFO_tuple_type& curr_tuple, void*& result_ptr, unsigned& num_result_elements)
{
  auto valid_result_found = false;
  auto query_field_idx = BCF_INFO_GET_QUERY_FIELD_IDX(curr_tuple);
  auto length_descriptor = m_query_config->get_length_descriptor_for_query_attribute_idx(query_field_idx);
  //Fields such as PL are skipped if the #alleles is above a certain threshold
  if(length_descriptor.is_length_genotype_dependent()
      && too_many_alt_alleles_for_genotype_length_fields(m_merged_alt_alleles.size()))
    return false;
  //Check if this is a field that was remapped - for remapped fields, we must use field objects from m_remapped_variant
  //else we should use field objects from the original variant
  auto& src_variant = (m_remapping_needed && (length_descriptor.is_length_allele_dependent()
        || query_field_idx == m_GT_query_idx))
    ? m_remapped_variant : variant;
  auto variant_type_enum = BCF_INFO_GET_VARIANT_FIELD_TYPE_ENUM(curr_tuple);
  //valid field handler
  assert(variant_type_enum < m_field_handlers.size() && m_field_handlers[variant_type_enum].get());
  auto num_valid_input_elements = 0u;
  switch(BCF_INFO_GET_VCF_FIELD_COMBINE_OPERATION(curr_tuple))
  {
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_SUM:
      valid_result_found = m_field_handlers[variant_type_enum]->get_valid_sum(src_variant, *m_query_config,
          query_field_idx, result_ptr, num_valid_input_elements);
      break;
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MEAN:
      valid_result_found = m_field_handlers[variant_type_enum]->get_valid_mean(src_variant, *m_query_config,
          query_field_idx, result_ptr, num_valid_input_elements);
      break;
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MEDIAN:
      valid_result_found = m_field_handlers[variant_type_enum]->get_valid_median(src_variant, *m_query_config,
          query_field_idx, result_ptr, num_valid_input_elements);
      break;
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_ELEMENT_WISE_SUM:
      valid_result_found = m_field_handlers[variant_type_enum]->compute_valid_element_wise_sum(src_variant, *m_query_config,
          query_field_idx, const_cast<const void**>(&result_ptr), num_result_elements);
      break;
    case VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_CONCATENATE:
      valid_result_found = m_field_handlers[variant_type_enum]->concatenate_field(src_variant, *m_query_config,
          query_field_idx, const_cast<const void**>(&result_ptr), num_result_elements);
      break;
    default:
      throw BroadCombinedGVCFException(std::string("Unknown VCF field combine operation ")
          +std::to_string(BCF_INFO_GET_VCF_FIELD_COMBINE_OPERATION(curr_tuple)));
      break;
  }
  return valid_result_found;
}

void BroadCombinedGVCFOperator::handle_INFO_fields(const Variant& variant)
{
  //interval variant, add END tag
  if(m_remapped_variant.get_column_end() > m_remapped_variant.get_column_begin())
  {
    int vcf_end_pos = m_remapped_variant.get_column_end() - m_curr_contig_begin_position + 1; //vcf END is 1 based
    bcf_update_info_int32(m_vcf_hdr, m_bcf_out, "END", &vcf_end_pos, 1);
    m_bcf_record_size += sizeof(int);
  }
  for(auto i=0u;i<m_INFO_fields_vec.size();++i)
  {
    auto& curr_tuple = m_INFO_fields_vec[i];
    //Just need a 4-byte value, the contents could be a float or int (determined by the templated median function)
    int32_t result = -1;
    void* result_ptr = reinterpret_cast<void*>(&result);
    //For element wise operations
    auto num_result_elements = 1u;
    auto valid_result_found = handle_VCF_field_combine_operation(variant, curr_tuple, result_ptr, num_result_elements);
    if(valid_result_found)
    {
      bcf_update_info(m_vcf_hdr, m_bcf_out, BCF_INFO_GET_VCF_FIELD_NAME(curr_tuple).c_str(), result_ptr, num_result_elements, BCF_INFO_GET_BCF_HT_TYPE(curr_tuple));
      m_bcf_record_size += num_result_elements*VariantFieldTypeUtil::size(BCF_INFO_GET_VARIANT_FIELD_TYPE_ENUM(curr_tuple));
    }
  }
}

void BroadCombinedGVCFOperator::handle_FORMAT_fields(const Variant& variant)
{
  //For weird DP field handling
  auto valid_DP_found = false;
  auto valid_MIN_DP_found = false; 
  m_MIN_DP_vector.resize(m_remapped_variant.get_num_calls()); //will do nothing after the first resize
  auto valid_DP_FORMAT_found = false; 
  m_DP_FORMAT_vector.resize(m_remapped_variant.get_num_calls());
  //Pointer to extended vector inside field handler object 
  const void* ptr = 0;
  uint64_t num_elements = 0ull;
  int* int_vec = 0;     //will point to same address as ptr, but as int*
  //Handle all fields - simply need to extend to the largest size
  for(auto i=0u;i<m_FORMAT_fields_vec.size();++i)
  {
    auto& curr_tuple = m_FORMAT_fields_vec[i];
    auto query_field_idx = BCF_FORMAT_GET_QUERY_FIELD_IDX(curr_tuple);
    auto length_descriptor = m_query_config->get_length_descriptor_for_query_attribute_idx(query_field_idx);
    //Fields such as PL are skipped if the #alleles is above a certain threshold
    if(length_descriptor.is_length_genotype_dependent()
        && too_many_alt_alleles_for_genotype_length_fields(m_merged_alt_alleles.size()))
      continue;
    auto known_field_enum = BCF_FORMAT_GET_KNOWN_FIELD_ENUM(curr_tuple);
    auto variant_type_enum = BCF_FORMAT_GET_VARIANT_FIELD_TYPE_ENUM(curr_tuple);
    auto is_char_type = (variant_type_enum == VARIANT_FIELD_CHAR);
    //valid field handler
    assert(variant_type_enum < m_field_handlers.size() && m_field_handlers[variant_type_enum].get());
    //Check if this is a field that was remapped - for remapped fields, we must use field objects from m_remapped_variant
    //else we should use field objects from the original variant
    auto& src_variant = (m_remapping_needed && (length_descriptor.is_length_allele_dependent()
          || query_field_idx == m_GT_query_idx))
        ? m_remapped_variant : variant;
    auto valid_field_found = m_field_handlers[variant_type_enum]->collect_and_extend_fields(src_variant, *m_query_config,
        query_field_idx, &ptr, num_elements,
        m_use_missing_values_not_vector_end && !is_char_type, m_use_missing_values_not_vector_end && is_char_type,
        known_field_enum == GVCF_GT_IDX);
    if(valid_field_found)
    {
      auto j=0u;
      auto do_insert = true;    //by default, insert into VCF record
      switch(known_field_enum)
      {
        case GVCF_GT_IDX: //GT field is a pita
          {
            int_vec = const_cast<int*>(reinterpret_cast<const int*>(ptr));
            auto num_elements_per_sample = num_elements/variant.get_num_calls();
            //GT of type 0/2 is stored as [0,0,2] in GenomicsDB if phased ploidy is used, else [0, 2]
            //GT of type 0|2 is stored as [0,1,2] in GenomicsDB if phased ploidy is used, else [0, 2]
            auto max_ploidy = length_descriptor.get_ploidy(num_elements_per_sample);
            uint64_t output_idx = 0ull;
            for(j=0ull;j<num_elements;j+=num_elements_per_sample)
              (*m_encode_GT_vector_function_ptr)(int_vec, j, num_elements_per_sample, output_idx);
            num_elements = max_ploidy*variant.get_num_calls();
            break;
          }
        case GVCF_GQ_IDX:
          do_insert = m_should_add_GQ_field;
          break;
        case GVCF_MIN_DP_IDX: //simply copy over min-dp values
          memcpy(&(m_MIN_DP_vector[0]), ptr, m_remapped_variant.get_num_calls()*sizeof(int));
          valid_MIN_DP_found = true;
          break;
        case GVCF_DP_FORMAT_IDX:
          memcpy(&(m_DP_FORMAT_vector[0]), ptr, m_remapped_variant.get_num_calls()*sizeof(int));
          valid_DP_FORMAT_found = true;
          do_insert = false; //Do not insert DP_FORMAT, wait till DP is read
          break;
        case GVCF_DP_IDX:
          int_vec = const_cast<int*>(reinterpret_cast<const int*>(ptr));
          valid_DP_found = true;
          do_insert = false; //Do not insert DP, handle DP related garbage at the end
          break;
        default:
          break;
      }
      if(do_insert)
      {
        bcf_update_format(m_vcf_hdr, m_bcf_out, BCF_FORMAT_GET_VCF_FIELD_NAME(curr_tuple).c_str(), ptr, num_elements,
            BCF_FORMAT_GET_BCF_HT_TYPE(curr_tuple));
        m_bcf_record_size += num_elements*VariantFieldTypeUtil::size(static_cast<VariantFieldTypeEnum>(variant_type_enum));
      }
    }
  }
  //Update DP fields
  if(valid_DP_found || valid_DP_FORMAT_found)
  {
    int sum_INFO_DP = 0;
    auto found_one_valid_DP_FORMAT = false;
    for(auto j=0ull;j<m_remapped_variant.get_num_calls();++j)
    {
      int dp_info_val = valid_DP_found ? int_vec[j] : get_bcf_missing_value<int>();
      int dp_format_val = valid_DP_FORMAT_found ? m_DP_FORMAT_vector[j] : get_bcf_missing_value<int>();
      assert(is_bcf_valid_value<int>(dp_format_val) || dp_format_val == get_bcf_missing_value<int>());
      assert(is_bcf_valid_value<int>(dp_info_val) || dp_info_val == get_bcf_missing_value<int>());
      //If DP_FORMAT is invalid, use dp_info value for DP_FORMAT
      //dp_format_val = is_bcf_valid_value<int>(dp_format_val) ? dp_format_val : dp_info_val;
      if(!is_bcf_valid_value<int>(dp_info_val))  //no valid DP info value found
      {
        //MIN_DP gets higher priority
        if(valid_MIN_DP_found && is_bcf_valid_value<int>(m_MIN_DP_vector[j]))
          dp_info_val = m_MIN_DP_vector[j];
        else //else DP_FORMAT
          dp_info_val = dp_format_val;
      }
      m_DP_FORMAT_vector[j] = dp_format_val;
      found_one_valid_DP_FORMAT = is_bcf_valid_value<int>(dp_format_val) || found_one_valid_DP_FORMAT;
      sum_INFO_DP += (is_bcf_valid_value<int>(dp_info_val) ? dp_info_val : 0);
    }
    if(found_one_valid_DP_FORMAT)
    {
      bcf_update_format_int32(m_vcf_hdr, m_bcf_out, "DP", &(m_DP_FORMAT_vector[0]), m_DP_FORMAT_vector.size()); //add DP FORMAT field
      m_bcf_record_size += m_DP_FORMAT_vector.size()*sizeof(int);
    }
    //If at least one valid DP value found from (DP or DP_FORMAT or MIN_DP), add DP to INFO
    if(sum_INFO_DP > 0 && !m_is_reference_block_only)
    {
      bcf_update_info_int32(m_vcf_hdr, m_bcf_out, "DP", &sum_INFO_DP, 1);
      m_bcf_record_size += sizeof(int);
    }
  }
}

//FIXME: totally naive implementation - too many reallocations etc
void BroadCombinedGVCFOperator::merge_ID_field(const Variant& variant, const unsigned query_idx)
{
  std::unordered_set<std::string> id_set;
  for(const auto& curr_call : variant)
  {
    auto& field_ptr = curr_call.get_field(query_idx);
    if(field_ptr.get() && field_ptr->is_valid())
    {
      auto* ptr = dynamic_cast<VariantFieldString*>(field_ptr.get());
      assert(ptr);
      const auto& curr_ID_value = ptr->get();
      auto last_begin_value = 0u;
      for(auto i=0u;i<curr_ID_value.length();++i)
        if(curr_ID_value[i] == ';')
        {
          id_set.insert(curr_ID_value.substr(last_begin_value, i-last_begin_value));
          last_begin_value = i+1u;
        }
      if(curr_ID_value.length() > last_begin_value)
        id_set.insert(curr_ID_value.substr(last_begin_value, curr_ID_value.length()-last_begin_value));
    }
  }
  m_ID_value.clear();
  for(const auto& str : id_set)
    m_ID_value += (str + ';');
  if(!m_ID_value.empty())
    m_ID_value.pop_back(); //delete last ';'
}

void BroadCombinedGVCFOperator::operate(Variant& variant, const VariantQueryConfig& query_config)
{
#ifdef DO_PROFILING
  m_bcf_t_creation_timer.start();
#endif
  //Handle spanning deletions - change ALT alleles in calls with deletions to *, <NON_REF>
  handle_deletions(variant, query_config);
  GA4GHOperator::operate(variant, query_config);
  //Moved to new contig
  if(static_cast<int64_t>(m_remapped_variant.get_column_begin()) >= m_next_contig_begin_position)
  {
    std::string contig_name;
    int64_t contig_position;
    auto status = m_vid_mapper->get_contig_location(m_remapped_variant.get_column_begin(), contig_name, contig_position);
    if(status)
    {
      int64_t contig_begin_position = m_remapped_variant.get_column_begin() - contig_position;
      if(contig_begin_position != m_next_contig_begin_position)
      {
        m_next_contig_name = std::move(contig_name);
        m_next_contig_begin_position = contig_begin_position;
      }
    }
    else
      throw BroadCombinedGVCFException("Unknown contig for position "+std::to_string(m_remapped_variant.get_column_begin()));
    switch_contig();
  }
  //clear out
  bcf_clear(m_bcf_out);
  m_bcf_record_size = 0ull;
  //If no valid FORMAT fields exist, this value is never set
  m_bcf_out->n_sample = bcf_hdr_nsamples(m_vcf_hdr);
  //position
  m_bcf_out->rid = m_curr_contig_hdr_idx;
  m_bcf_out->pos = m_remapped_variant.get_column_begin() - m_curr_contig_begin_position;
  //ID field
  if(m_query_config->is_defined_query_idx_for_known_field_enum(GVCF_ID_IDX))
  {
    auto ID_query_idx = m_query_config->get_query_idx_for_known_field_enum(GVCF_ID_IDX);
    merge_ID_field(variant, ID_query_idx);
    if(!m_ID_value.empty())
      bcf_update_id(m_vcf_hdr, m_bcf_out, m_ID_value.c_str());
  }
  m_bcf_record_size += m_ID_value.length();
  //GATK combined GVCF does not care about QUAL value
  m_bcf_out->qual = get_bcf_missing_value<float>();
  if(BCF_INFO_GET_VCF_FIELD_COMBINE_OPERATION(m_vcf_qual_tuple) != VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_UNKNOWN_OPERATION)
  {
    unsigned num_result_elements = 1u;
    auto qual_result = 1.0f;
    void* result_ptr = reinterpret_cast<void*>(&qual_result);
    auto valid_result_found = handle_VCF_field_combine_operation(variant, m_vcf_qual_tuple, result_ptr, num_result_elements);
    if(valid_result_found)
      m_bcf_out->qual = qual_result;
  }
  m_bcf_record_size += 3*sizeof(int);
  //Update alleles
  auto& ref_allele = dynamic_cast<VariantFieldString*>(m_remapped_variant.get_common_field(0u).get())->get();
  if(ref_allele.length() == 1u && ref_allele[0] == 'N')
  {
    ref_allele[0] = m_vcf_adapter->get_reference_base_at_position(m_curr_contig_name.c_str(), m_bcf_out->pos);
    if(BroadCombinedGVCFOperator::m_legal_bases.find(ref_allele[0]) == BroadCombinedGVCFOperator::m_legal_bases.end())
      ref_allele[0] = 'N';
  }
  const auto& alt_alleles = dynamic_cast<VariantFieldALTData*>(m_remapped_variant.get_common_field(1u).get())->get();
  auto total_num_merged_alleles = alt_alleles.size() + 1u;      //+1 for REF
  if(total_num_merged_alleles > m_alleles_pointer_buffer.size())
    m_alleles_pointer_buffer.resize(total_num_merged_alleles);
  //REF
  m_alleles_pointer_buffer[0] = ref_allele.c_str();
  m_bcf_record_size += ref_allele.length()*sizeof(char);
  //ALT
  for(auto i=1u;i<total_num_merged_alleles;++i)
  {
    m_alleles_pointer_buffer[i] = alt_alleles[i-1u].c_str();
    m_bcf_record_size += alt_alleles[i-1u].length()*sizeof(char);
  }
  bcf_update_alleles(m_vcf_hdr, m_bcf_out, &(m_alleles_pointer_buffer[0]), total_num_merged_alleles);
  //FILTER fields
  if(m_vcf_adapter->produce_FILTER_field() &&
      m_query_config->is_defined_query_idx_for_known_field_enum(GVCF_FILTER_IDX))
  {
    auto FILTER_query_idx = m_query_config->get_query_idx_for_known_field_enum(GVCF_FILTER_IDX);
    //Remove duplicates across samples
    std::unordered_set<int> filter_idx_set;
    for(const auto& call : variant)
    {
      auto curr_FILTER_field = call.get_field<VariantFieldPrimitiveVectorData<int>>(FILTER_query_idx);
      if(curr_FILTER_field && curr_FILTER_field->is_valid())
      {
        auto& curr_FILTER_idx_vec = curr_FILTER_field->get();
        filter_idx_set.insert(curr_FILTER_idx_vec.begin(), curr_FILTER_idx_vec.end());
      }
    }
    if(filter_idx_set.size())
    {
      m_FILTER_idx_vec.resize(filter_idx_set.size());
      auto idx = 0u;
      for(auto global_field_idx : filter_idx_set)
      {
        assert(static_cast<size_t>(global_field_idx) < m_global_field_idx_to_hdr_idx.size());
        auto hdr_field_idx = m_global_field_idx_to_hdr_idx[global_field_idx];
        if(hdr_field_idx >= 0)
          m_FILTER_idx_vec[idx++] = hdr_field_idx;
      }
      bcf_update_filter(m_vcf_hdr, m_bcf_out, &(m_FILTER_idx_vec[0]), m_FILTER_idx_vec.size());
    }
  }
  //Flag that determines when to add GQ field - only when <NON_REF> is the only alternate allele
  //m_should_add_GQ_field = (m_NON_REF_exists && alt_alleles.size() == 1u);
  m_should_add_GQ_field = true; //always added in new version of CombineGVCFs
  //INFO fields
  handle_INFO_fields(variant);
  //FORMAT fields
  handle_FORMAT_fields(variant);
#ifdef DO_PROFILING
  m_bcf_t_creation_timer.stop();
#endif
  m_vcf_adapter->handoff_output_bcf_line(m_bcf_out, m_bcf_record_size);
}

void BroadCombinedGVCFOperator::switch_contig()
{
  m_curr_contig_name = std::move(m_next_contig_name);
  m_curr_contig_begin_position = m_next_contig_begin_position;
  m_curr_contig_hdr_idx = bcf_hdr_id2int(m_vcf_hdr, BCF_DT_CTG, m_curr_contig_name.c_str());
  m_vid_mapper->get_next_contig_location(m_next_contig_begin_position, m_next_contig_name, m_next_contig_begin_position);
}

//Modifies original Variant object
void BroadCombinedGVCFOperator::handle_deletions(Variant& variant, const VariantQueryConfig& query_config)
{
  m_reduced_alleles_LUT.resize_luts_if_needed(variant.get_num_calls(), 10u);    //will not have more than 3 alleles anyway
  m_reduced_alleles_LUT.reset_luts();
  auto GT_length_descriptor = m_query_config->get_length_descriptor_for_query_attribute_idx(m_GT_query_idx);
  for(auto iter=variant.begin(), e=variant.end();iter != e;++iter)
  {
    auto& curr_call = *iter;
    auto curr_call_idx_in_variant = iter.get_call_idx_in_variant();
    //Deletion and not handled as spanning deletion 
    //So replace the ALT with *,<NON_REF> and REF with "N"
    //Remap PL, AD fields
    if(curr_call.contains_deletion() && variant.get_column_begin() > curr_call.get_column_begin())
    {
      auto& ref_allele = get_known_field<VariantFieldString, true>(curr_call, query_config, GVCF_REF_IDX)->get();
      auto& alt_alleles = get_known_field<VariantFieldALTData, true>(curr_call, query_config, GVCF_ALT_IDX)->get();
      assert(alt_alleles.size() > 0u);
      //Already handled as a spanning deletion, nothing to do
      if(alt_alleles[0u] == g_vcf_SPANNING_DELETION)
        continue;
      //Reduced allele list will be REF="N", ALT="*, <NON_REF>"
      m_reduced_alleles_LUT.add_input_merged_idx_pair(curr_call_idx_in_variant, 0, 0);  //REF-REF mapping
      //For each deletion allele, find the PL value corresponding to the genotype in which all the elements
      //of the genotype are equal to the deletion allele. The deletion allele with the lowest PL value is
      //mapped to "*" allele
      //GT field - for ploidy
      auto ploidy = 0u;
      auto* original_GT_field_ptr = (m_GT_query_idx != UNDEFINED_ATTRIBUTE_IDX_VALUE)
          ? curr_call.get_field<VariantFieldPrimitiveVectorData<int>>(m_GT_query_idx) : 0;
      if(original_GT_field_ptr && original_GT_field_ptr->is_valid())
        ploidy = GT_length_descriptor.get_ploidy(original_GT_field_ptr->get().size());
      else
        original_GT_field_ptr = 0;
      auto lowest_deletion_allele_idx = -1;
      int lowest_PL_value = INT_MAX;
      auto PL_field_ptr = get_known_field_if_queried<VariantFieldPrimitiveVectorData<int>, true>(curr_call, query_config, GVCF_PL_IDX);
      auto has_NON_REF = false;
      //PL field exists
      if(PL_field_ptr && PL_field_ptr->is_valid())
      {
        auto& PL_vector = PL_field_ptr->get();
        m_spanning_deletion_current_genotype.resize(ploidy);
        for(auto i=0u;i<alt_alleles.size();++i)
        {
          auto allele_idx = i+1;  //+1 for REF
          if(VariantUtils::is_deletion(ref_allele, alt_alleles[i]))
          {
            if(lowest_deletion_allele_idx < 0) //uninitialized
              lowest_deletion_allele_idx = allele_idx;
            //Genotype with all elements set to the deletion allele
            m_spanning_deletion_current_genotype.assign(ploidy, allele_idx);
            auto gt_idx = VariantOperations::get_genotype_index(m_spanning_deletion_current_genotype, true);
            //Truncations - dropped values etc
            if(gt_idx < PL_vector.size() && PL_vector[gt_idx] < lowest_PL_value)
            {
              lowest_PL_value = PL_vector[gt_idx];
              lowest_deletion_allele_idx = allele_idx;
            }
          }
          else
            if(IS_NON_REF_ALLELE(alt_alleles[i]))
            {
              m_reduced_alleles_LUT.add_input_merged_idx_pair(curr_call_idx_in_variant, allele_idx, 2);
              has_NON_REF = true;
            }
        }
      }
      else      //PL field is not queried, simply use the first ALT allele
        lowest_deletion_allele_idx = 1;
      assert(lowest_deletion_allele_idx >= 1);    //should be an ALT allele
      //first ALT allele in reduced list is *
      m_reduced_alleles_LUT.add_input_merged_idx_pair(curr_call_idx_in_variant, lowest_deletion_allele_idx, 1); 
      if(has_NON_REF)
      {
        alt_alleles.resize(2u);
        alt_alleles[1u] = TILEDB_NON_REF_VARIANT_REPRESENTATION;
      }
      else
        alt_alleles.resize(1u); //only spanning deletion
      ref_allele = "N"; //set to unknown REF for now
      alt_alleles[0u] = g_vcf_SPANNING_DELETION;
      unsigned num_reduced_alleles = alt_alleles.size() + 1u;   //+1 for REF
      //GT field
      if(original_GT_field_ptr)
      {
        auto& input_GT =
          original_GT_field_ptr->get();
        m_spanning_deletion_remapped_GT.resize(input_GT.size());
        VariantOperations::remap_GT_field(input_GT, m_spanning_deletion_remapped_GT, m_reduced_alleles_LUT, curr_call_idx_in_variant,
            num_reduced_alleles, has_NON_REF, GT_length_descriptor);
        //Copy back
        memcpy(&(input_GT[0]), &(m_spanning_deletion_remapped_GT[0]), input_GT.size()*sizeof(int));
      }
      //Remap fields that need to be remapped
      for(auto i=0u;i<m_remapped_fields_query_idxs.size();++i)
      {
        auto query_field_idx = m_remapped_fields_query_idxs[i];
        auto length_descriptor = query_config.get_length_descriptor_for_query_attribute_idx(query_field_idx);
        //field whose length is dependent on #alleles
        assert(length_descriptor.is_length_allele_dependent());
        unsigned num_reduced_elements =
          length_descriptor.get_num_elements(num_reduced_alleles-1u, ploidy, 0u);     //#alt alleles
        //Remapper for variant
        RemappedVariant remapper_variant(variant, query_field_idx); 
        auto& curr_field = curr_call.get_field(query_field_idx);
        if(curr_field.get() && curr_field->is_valid())      //Not null
        {
          //Copy field to pass to remap function 
          assert(i < m_spanning_deletions_remapped_fields.size());
          copy_field(m_spanning_deletions_remapped_fields[i], curr_field);
          curr_field->resize(num_reduced_elements);
          //Get handler for current type
          auto& handler = get_handler_for_type(query_config.get_element_type(query_field_idx));
          assert(handler.get());
          //Call remap function
          handler->remap_vector_data(
              m_spanning_deletions_remapped_fields[i], curr_call_idx_in_variant,
              m_reduced_alleles_LUT, num_reduced_alleles, has_NON_REF, ploidy,
              length_descriptor, num_reduced_elements, remapper_variant);
        }
      }
      //Invalidate INFO fields
      for(const auto& tuple : m_INFO_fields_vec)
      {
        auto query_idx = BCF_INFO_GET_QUERY_FIELD_IDX(tuple);
        auto& field = curr_call.get_field(query_idx);
        if(field.get())
          field->set_valid(false);
      }
    }
  }
}

#endif //ifdef HTSDIR
