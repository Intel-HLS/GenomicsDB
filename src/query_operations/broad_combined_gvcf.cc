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

#ifdef HTSDIR

#include "broad_combined_gvcf.h"

//INFO fields
#define MAKE_BCF_INFO_TUPLE(enum_idx, variant_type_enum, bcf_type) \
  INFO_tuple_type(enum_idx, variant_type_enum, bcf_type)
#define BCF_INFO_GET_KNOWN_FIELD_ENUM(X) (std::get<0>(X))
#define BCF_INFO_GET_VARIANT_FIELD_TYPE_ENUM(X) (std::get<1>(X))
#define BCF_INFO_GET_BCF_HT_TYPE(X) (std::get<2>(X))
//FORMAT fields
#define MAKE_BCF_FORMAT_TUPLE(enum_idx, variant_type_enum, bcf_type) \
  FORMAT_tuple_type(enum_idx, variant_type_enum, bcf_type)
#define BCF_FORMAT_GET_KNOWN_FIELD_ENUM(X) (std::get<0>(X))
#define BCF_FORMAT_GET_QUERY_FIELD_IDX(X) (std::get<0>(X))
#define BCF_FORMAT_GET_VARIANT_FIELD_TYPE_ENUM(X) (std::get<1>(X))
#define BCF_FORMAT_GET_BCF_HT_TYPE(X) (std::get<2>(X))

//Static member
const std::unordered_set<char> BroadCombinedGVCFOperator::m_legal_bases({'A', 'T', 'G', 'C'});

BroadCombinedGVCFOperator::BroadCombinedGVCFOperator(VCFAdapter& vcf_adapter, const VidMapper& id_mapper,
    const VariantQueryConfig& query_config,
    const unsigned max_diploid_alt_alleles_that_can_be_genotyped, const bool use_missing_values_only_not_vector_end)
: GA4GHOperator(query_config, max_diploid_alt_alleles_that_can_be_genotyped)
{
  clear();
  if(!id_mapper.is_initialized())
    throw BroadCombinedGVCFException("Id mapper is not initialized");
  m_query_config = &query_config;
  //Initialize VCF structs
  m_vcf_adapter = &vcf_adapter;
  m_vid_mapper = &id_mapper;
  m_use_missing_values_not_vector_end = use_missing_values_only_not_vector_end;
  m_vcf_hdr = vcf_adapter.get_vcf_header();
  m_bcf_out = bcf_init();
  //vector of char*, to avoid frequent reallocs()
  m_alleles_pointer_buffer.resize(100u);
  //INFO fields
  m_INFO_fields_vec = std::move(std::vector<INFO_tuple_type>
      {
      MAKE_BCF_INFO_TUPLE(GVCF_BASEQRANKSUM_IDX, VARIANT_FIELD_FLOAT, BCF_HT_REAL),
      MAKE_BCF_INFO_TUPLE(GVCF_CLIPPINGRANKSUM_IDX, VARIANT_FIELD_FLOAT, BCF_HT_REAL),
      MAKE_BCF_INFO_TUPLE(GVCF_MQRANKSUM_IDX, VARIANT_FIELD_FLOAT, BCF_HT_REAL),
      MAKE_BCF_INFO_TUPLE(GVCF_READPOSRANKSUM_IDX, VARIANT_FIELD_FLOAT, BCF_HT_REAL),
      MAKE_BCF_INFO_TUPLE(GVCF_MQ_IDX, VARIANT_FIELD_FLOAT, BCF_HT_REAL),
      MAKE_BCF_INFO_TUPLE(GVCF_RAW_MQ_IDX, VARIANT_FIELD_FLOAT, BCF_HT_REAL),
      MAKE_BCF_INFO_TUPLE(GVCF_MQ0_IDX, VARIANT_FIELD_INT, BCF_HT_INT)
      });
  m_FORMAT_fields_vec = std::move(std::vector<FORMAT_tuple_type>
      {
      MAKE_BCF_FORMAT_TUPLE(GVCF_GT_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_GQ_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_SB_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_AD_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_PL_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_PGT_IDX, VARIANT_FIELD_CHAR, BCF_HT_STR),
      MAKE_BCF_FORMAT_TUPLE(GVCF_PID_IDX, VARIANT_FIELD_CHAR, BCF_HT_STR),
      MAKE_BCF_FORMAT_TUPLE(GVCF_MIN_DP_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_DP_FORMAT_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_DP_IDX, VARIANT_FIELD_INT, BCF_HT_INT)  //always last field, read DP not DP_FORMAT field
      });
  //Last field should always be DP
  if(BCF_FORMAT_GET_KNOWN_FIELD_ENUM(m_FORMAT_fields_vec[m_FORMAT_fields_vec.size()-1u]) != GVCF_DP_IDX)
    throw BroadCombinedGVCFException("Last queried FORMAT field should be DP, instead it is "
        +g_known_variant_field_names[BCF_FORMAT_GET_KNOWN_FIELD_ENUM(m_FORMAT_fields_vec[m_FORMAT_fields_vec.size()-1u])]);
  //Discard fields not part of the query
  auto last_valid_idx = 0u;
  for(auto i=0u;i<m_INFO_fields_vec.size();++i)
  {
    auto& tuple = m_INFO_fields_vec[i];
    if(query_config.is_defined_query_idx_for_known_field_enum((BCF_INFO_GET_KNOWN_FIELD_ENUM(tuple))))
    {
      m_INFO_fields_vec[last_valid_idx++] = tuple;
      VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, &id_mapper, KnownFieldInfo::get_known_field_name_for_enum(BCF_INFO_GET_KNOWN_FIELD_ENUM(tuple)),
          BCF_HL_INFO);
    }
  }
  //Add DP field to header
  VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, &id_mapper, "DP", BCF_HL_INFO);
  m_INFO_fields_vec.resize(last_valid_idx);
  //Same for FORMAT
  last_valid_idx = 0u;
  std::unordered_set<unsigned> handled_format_fields_query_idxs;
  for(auto i=0u;i<m_FORMAT_fields_vec.size();++i)
  {
    auto& tuple = m_FORMAT_fields_vec[i];
    auto known_field_enum = BCF_FORMAT_GET_KNOWN_FIELD_ENUM(tuple);
    if(query_config.is_defined_query_idx_for_known_field_enum(known_field_enum))
    {
      m_FORMAT_fields_vec[last_valid_idx++] = tuple;
      auto field_info_ptr = id_mapper.get_field_info(KnownFieldInfo::get_known_field_name_for_enum(known_field_enum));
      assert(field_info_ptr);
      VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, &id_mapper, field_info_ptr->m_name, BCF_HL_FMT);
      auto query_field_idx = query_config.get_query_idx_for_known_field_enum(known_field_enum);
      handled_format_fields_query_idxs.insert(query_field_idx);
    }
  }
  m_FORMAT_fields_vec.resize(last_valid_idx);
  //Add format fields which are not already known by GenomicsDB
  for(auto i=0u;i<query_config.get_num_queried_attributes();++i)
  {
    auto* field_info = m_vid_mapper->get_field_info(query_config.get_query_attribute_name(i));
    if(field_info &&
        field_info->m_is_vcf_FORMAT_field && handled_format_fields_query_idxs.find(i) == handled_format_fields_query_idxs.end())
    {
      m_unknown_FORMAT_fields_vec.emplace_back(MAKE_BCF_FORMAT_TUPLE(i,
            VariantFieldTypeUtil::get_variant_field_type_enum_for_variant_field_type(field_info->m_type_index),
            VariantFieldTypeUtil::get_vcf_field_type_enum_for_variant_field_type(field_info->m_type_index)));
      VCFAdapter::add_field_to_hdr_if_missing(m_vcf_hdr, &id_mapper, field_info->m_name, BCF_HL_FMT);
    }
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
    assert(status);
    bcf_hdr_add_sample(m_vcf_hdr, callset_name.c_str());
  }
  bcf_hdr_sync(m_vcf_hdr);
  m_vcf_adapter->print_header();
  //vector of field pointers used for handling remapped fields when dealing with spanning deletions
  //Individual pointers will be allocated later
  m_spanning_deletions_remapped_fields.resize(m_remapped_fields_query_idxs.size());
}

void BroadCombinedGVCFOperator::clear()
{
  m_curr_contig_name.clear();
  m_next_contig_name.clear();
  m_alleles_pointer_buffer.clear();
  m_INFO_fields_vec.clear();
  m_FORMAT_fields_vec.clear();
  m_unknown_FORMAT_fields_vec.clear();
  m_MIN_DP_vector.clear();
  m_DP_FORMAT_vector.clear();
  m_spanning_deletions_remapped_fields.clear();
}

void BroadCombinedGVCFOperator::handle_INFO_fields(const Variant& variant)
{
  //interval variant, add END tag
  if(m_remapped_variant.get_column_end() > m_remapped_variant.get_column_begin())
  {
    int vcf_end_pos = m_remapped_variant.get_column_end() - m_curr_contig_begin_position + 1; //vcf END is 1 based
    bcf_update_info_int32(m_vcf_hdr, m_bcf_out, "END", &vcf_end_pos, 1);
  }
  //Compute median for all INFO fields except RAW_MQ
  for(auto i=0u;i<m_INFO_fields_vec.size();++i)
  {
    auto& curr_tuple = m_INFO_fields_vec[i];
    auto known_field_enum = BCF_INFO_GET_KNOWN_FIELD_ENUM(curr_tuple);
    assert(known_field_enum < g_known_variant_field_names.size());
    auto variant_type_enum = BCF_INFO_GET_VARIANT_FIELD_TYPE_ENUM(curr_tuple);
    //valid field handler
    assert(variant_type_enum < m_field_handlers.size() && m_field_handlers[variant_type_enum].get());
    if(known_field_enum != GVCF_RAW_MQ_IDX)
    {
      //Just need a 4-byte value, the contents could be a float or int (determined by the templated median function)
      int32_t median;
      auto valid_median_found = m_field_handlers[variant_type_enum]->get_valid_median(variant, *m_query_config,
          m_query_config->get_query_idx_for_known_field_enum(known_field_enum), reinterpret_cast<void*>(&median));
      if(valid_median_found)
        bcf_update_info(m_vcf_hdr, m_bcf_out, g_known_variant_field_names[known_field_enum].c_str(), &median, 1, BCF_INFO_GET_BCF_HT_TYPE(curr_tuple));
    }
    else
    {
      //RAW_MQ - was element wise sum initially
      //auto num_elements = KnownFieldInfo::get_num_elements_for_known_field_enum(GVCF_RAW_MQ_IDX, 0u, 0u);
      //auto result_vector = std::move(std::vector<int>(num_elements, 0u));
      //auto valid_sum_found = m_field_handlers[variant_type_enum]->compute_valid_element_wise_sum(variant, *m_query_config,
      //m_query_config->get_query_idx_for_known_field_enum(known_field_enum), reinterpret_cast<void*>(&(result_vector[0])), num_elements);
      //Just need a 4-byte value, the contents could be a float or int (determined by the templated sum function)
      int32_t sum;
      auto valid_sum_found = m_field_handlers[variant_type_enum]->get_valid_sum(variant, *m_query_config,
          m_query_config->get_query_idx_for_known_field_enum(known_field_enum), reinterpret_cast<void*>(&sum));
      if(valid_sum_found)
        bcf_update_info(m_vcf_hdr, m_bcf_out, g_known_variant_field_names[known_field_enum].c_str(), &sum, 1, BCF_INFO_GET_BCF_HT_TYPE(curr_tuple));
      //bcf_update_info(m_vcf_hdr, m_bcf_out, g_known_variant_field_names[known_field_enum].c_str(), &(result_vector[0]), num_elements,
      //BCF_INFO_GET_BCF_HT_TYPE(curr_tuple));
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
  auto num_elements = 0u;
  int* int_vec = 0;     //will point to same address as ptr, but as int*
  //Handle all fields - simply need to extend to the largest size
  for(auto i=0u;i<m_FORMAT_fields_vec.size();++i)
  {
    auto& curr_tuple = m_FORMAT_fields_vec[i];
    auto known_field_enum = BCF_FORMAT_GET_KNOWN_FIELD_ENUM(curr_tuple);
    assert(known_field_enum < g_known_variant_field_names.size());
    if(KnownFieldInfo::is_length_genotype_dependent(known_field_enum) && too_many_alt_alleles_for_genotype_length_fields(m_merged_alt_alleles.size()))
      continue;
    auto variant_type_enum = BCF_FORMAT_GET_VARIANT_FIELD_TYPE_ENUM(curr_tuple);
    auto is_char_type = (variant_type_enum == VARIANT_FIELD_CHAR);
    //valid field handler
    assert(variant_type_enum < m_field_handlers.size() && m_field_handlers[variant_type_enum].get());
    //Check if this is a field that was remapped - for remapped fields, we must use field objects from m_remapped_variant
    //else we should use field objects from the original variant
    auto query_field_idx = m_query_config->get_query_idx_for_known_field_enum(known_field_enum);
    auto& src_variant = (m_remapping_needed && KnownFieldInfo::is_length_allele_dependent(known_field_enum)) ? m_remapped_variant : variant;
    auto valid_field_found = m_field_handlers[variant_type_enum]->collect_and_extend_fields(src_variant, *m_query_config,
        query_field_idx, &ptr, num_elements,
        m_use_missing_values_not_vector_end && !is_char_type, m_use_missing_values_not_vector_end && is_char_type);
    if(valid_field_found)
    {
      auto j=0u;
      auto do_insert = true;    //by default, insert into VCF record
      switch(known_field_enum)
      {
        case GVCF_GT_IDX: //GT field is a pita
          int_vec = const_cast<int*>(reinterpret_cast<const int*>(ptr));
          //CombineGVCF sets GT field to missing
          for(j=0u;j<num_elements;++j)
            int_vec[j] = bcf_gt_missing;
          //int_vec[j] = (int_vec[j] == get_bcf_missing_value<int>()) ? bcf_gt_missing : bcf_gt_unphased(int_vec[j]);
          break;
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
        bcf_update_format(m_vcf_hdr, m_bcf_out, g_known_variant_field_names[known_field_enum].c_str(), ptr, num_elements,
            BCF_FORMAT_GET_BCF_HT_TYPE(curr_tuple));
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
      bcf_update_format_int32(m_vcf_hdr, m_bcf_out, "DP", &(m_DP_FORMAT_vector[0]), m_DP_FORMAT_vector.size()); //add DP FORMAT field
    //If at least one valid DP value found from (DP or DP_FORMAT or MIN_DP), add DP to INFO
    if(sum_INFO_DP > 0 && !m_is_reference_block_only)
      bcf_update_info_int32(m_vcf_hdr, m_bcf_out, "DP", &sum_INFO_DP, 1);
  }
  //Handle fields which GenomicsDB does not know about 
  for(auto i=0u;i<m_unknown_FORMAT_fields_vec.size();++i)
  {
    auto& curr_tuple = m_unknown_FORMAT_fields_vec[i];
    auto query_field_idx = BCF_FORMAT_GET_QUERY_FIELD_IDX(curr_tuple);
    auto variant_type_enum = BCF_FORMAT_GET_VARIANT_FIELD_TYPE_ENUM(curr_tuple);
    auto is_char_type = (variant_type_enum == VARIANT_FIELD_CHAR);
    //valid field handler
    assert(variant_type_enum < m_field_handlers.size() && m_field_handlers[variant_type_enum].get());
    auto valid_field_found = m_field_handlers[variant_type_enum]->collect_and_extend_fields(variant, *m_query_config,
        query_field_idx, &ptr, num_elements,
        m_use_missing_values_not_vector_end && !is_char_type, m_use_missing_values_not_vector_end && is_char_type);
    if(valid_field_found)
      bcf_update_format(m_vcf_hdr, m_bcf_out, m_query_config->get_query_attribute_name(query_field_idx).c_str(), ptr, num_elements,
          BCF_FORMAT_GET_BCF_HT_TYPE(curr_tuple));
  }
}

void BroadCombinedGVCFOperator::operate(Variant& variant, const VariantQueryConfig& query_config)
{
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
  //If no valid FORMAT fields exist, this value is never set
  m_bcf_out->n_sample = bcf_hdr_nsamples(m_vcf_hdr);
  //position
  m_bcf_out->rid = m_curr_contig_hdr_idx;
  m_bcf_out->pos = m_remapped_variant.get_column_begin() - m_curr_contig_begin_position;
  //GATK combined GVCF does not care about QUAL value
  m_bcf_out->qual = get_bcf_missing_value<float>();
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
  //ALT
  for(auto i=1u;i<total_num_merged_alleles;++i)
    m_alleles_pointer_buffer[i] = alt_alleles[i-1u].c_str();
  bcf_update_alleles(m_vcf_hdr, m_bcf_out, &(m_alleles_pointer_buffer[0]), total_num_merged_alleles);
  //Flag that determines when to add GQ field - only when <NON_REF> is the only alternate allele
  //m_should_add_GQ_field = (m_NON_REF_exists && alt_alleles.size() == 1u);
  m_should_add_GQ_field = true; //always added in new version of CombineGVCFs
  //INFO fields
  handle_INFO_fields(variant);
  //FORMAT fields
  handle_FORMAT_fields(variant);
  m_vcf_adapter->handoff_output_bcf_line(m_bcf_out);
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
      //Need to find deletion allele with lowest PL value - this deletion allele is mapped to "*" allele
      auto lowest_deletion_allele_idx = -1;
      int lowest_PL_value = INT_MAX;
      auto PL_field_ptr = get_known_field_if_queried<VariantFieldPrimitiveVectorData<int>, true>(curr_call, query_config, GVCF_PL_IDX);
      auto has_NON_REF = false;
      //PL field exists
      if(PL_field_ptr)
      {
        auto& PL_vector = PL_field_ptr->get();
        for(auto i=0u;i<alt_alleles.size();++i)
        {
          auto allele_idx = i+1;  //+1 for REF
          if(VariantUtils::is_deletion(ref_allele, alt_alleles[i]))
          {
            unsigned gt_idx = bcf_alleles2gt(allele_idx, allele_idx);
            assert(gt_idx < PL_vector.size());
            if(PL_vector[gt_idx] < lowest_PL_value)
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
      //Remap fields that need to be remapped
      for(auto i=0u;i<m_remapped_fields_query_idxs.size();++i)
      {
        auto query_field_idx = m_remapped_fields_query_idxs[i];
        //is known field?
        assert(query_config.is_defined_known_field_enum_for_query_idx(query_field_idx));
        const auto* info_ptr = query_config.get_info_for_query_idx(query_field_idx);
        //known field whose length is dependent on #alleles
        assert(info_ptr && info_ptr->is_length_allele_dependent());
        unsigned num_reduced_elements = info_ptr->get_num_elements_for_known_field_enum(num_reduced_alleles-1u, 0u);     //#alt alleles
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
          auto& handler = get_handler_for_type(curr_field->get_element_type());
          assert(handler.get());
          //Call remap function
          handler->remap_vector_data(
              m_spanning_deletions_remapped_fields[i], curr_call_idx_in_variant,
              m_reduced_alleles_LUT, num_reduced_alleles, has_NON_REF,
              info_ptr->get_length_descriptor(), num_reduced_elements, remapper_variant);
        }
      }
      //Broad's CombineGVCF ignores GT field anyway - set missing
      if(m_GT_query_idx != UNDEFINED_ATTRIBUTE_IDX_VALUE)
      {
        auto& GT_vector = get_known_field<VariantFieldPrimitiveVectorData<int>, true>(curr_call, query_config, GVCF_GT_IDX)->get();
        for(auto i=0u;i<GT_vector.size();++i)
          GT_vector[i] = get_tiledb_null_value<int>(); 
      }
      //Invalidate INFO fields
      for(const auto& tuple : m_INFO_fields_vec)
      {
        assert(query_config.is_defined_query_idx_for_known_field_enum(BCF_INFO_GET_KNOWN_FIELD_ENUM(tuple)));
        auto query_idx = query_config.get_query_idx_for_known_field_enum(BCF_INFO_GET_KNOWN_FIELD_ENUM(tuple));
        auto& field = curr_call.get_field(query_idx);
        if(field.get())
          field->set_valid(false);
      }
    }
  }
}

#endif //ifdef HTSDIR
