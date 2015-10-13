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
#define BCF_FORMAT_GET_VARIANT_FIELD_TYPE_ENUM(X) (std::get<1>(X))
#define BCF_FORMAT_GET_BCF_HT_TYPE(X) (std::get<2>(X))


BroadCombinedGVCFOperator::BroadCombinedGVCFOperator(VCFAdapter& vcf_adapter, const VariantQueryConfig& query_config) 
: GA4GHOperator()
{
  clear();
  m_query_config = &query_config;
  //Initialize VCF structs
  m_vcf_adapter = &vcf_adapter;
  m_vcf_hdr = vcf_adapter.get_vcf_header();
  m_bcf_out = bcf_init();
  //Get contig info for position 0, store curr contig in next_contig and call switch_contig function to do all the setup
  auto curr_contig_flag = m_vcf_adapter->get_contig_location(0, m_next_contig_name, m_next_contig_begin_position);
  assert(curr_contig_flag);
  switch_contig();
  //GATK combined GVCF does not care about QUAL value
  m_bcf_out->qual = get_bcf_missing_value<float>();
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
      MAKE_BCF_INFO_TUPLE(GVCF_MQ0_IDX, VARIANT_FIELD_INT, BCF_HT_INT)
      });
  m_FORMAT_fields_vec = std::move(std::vector<FORMAT_tuple_type>
      {
      MAKE_BCF_FORMAT_TUPLE(GVCF_GT_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_GQ_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_SB_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_AD_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_PL_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_MIN_DP_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_DP_FORMAT_IDX, VARIANT_FIELD_INT, BCF_HT_INT),
      MAKE_BCF_FORMAT_TUPLE(GVCF_DP_IDX, VARIANT_FIELD_INT, BCF_HT_INT)  //always last field, read DP not DP_FORMAT field
      });
  //Last field should always be DP
  if(BCF_FORMAT_GET_KNOWN_FIELD_ENUM(m_FORMAT_fields_vec[m_FORMAT_fields_vec.size()-1u]) != GVCF_DP_IDX)
    throw BroadCombinedGVCFException("Last queried FORMAT field should be DP, instead it is "
        +g_known_variant_field_names[BCF_FORMAT_GET_KNOWN_FIELD_ENUM(m_FORMAT_fields_vec[m_FORMAT_fields_vec.size()-1u])]);
  //Sanity checks - all required fields must be queried
  for(auto& tuple : m_INFO_fields_vec)
    if(!query_config.is_defined_query_idx_for_known_field_enum((BCF_INFO_GET_KNOWN_FIELD_ENUM(tuple))))
      throw BroadCombinedGVCFException("Field "+g_known_variant_field_names[BCF_INFO_GET_KNOWN_FIELD_ENUM(tuple)]+" not specified as part of query");
  for(auto& tuple : m_FORMAT_fields_vec)
    if(!query_config.is_defined_query_idx_for_known_field_enum((BCF_FORMAT_GET_KNOWN_FIELD_ENUM(tuple))))
      throw BroadCombinedGVCFException("Field "+g_known_variant_field_names[BCF_FORMAT_GET_KNOWN_FIELD_ENUM(tuple)]+" not specified as part of query");
  //Add samples to template header
  for(auto i=0ull;i<query_config.get_num_rows_to_query();++i)
  {
    auto row_idx = query_config.get_array_row_idx_for_query_row_idx(i);
    bcf_hdr_add_sample(m_vcf_hdr, m_vcf_adapter->get_sample_name_for_idx(row_idx));
  }
  bcf_hdr_sync(m_vcf_hdr);
  m_vcf_adapter->print_header();
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
}

void BroadCombinedGVCFOperator::handle_INFO_fields()
{
  auto& copy_variant = m_variants[0];
  //interval variant, add END tag
  if(copy_variant.get_column_end() > copy_variant.get_column_begin())
  {
    int vcf_end_pos = copy_variant.get_column_end() - m_curr_contig_begin_position + 1; //vcf END is 1 based
    bcf_update_info_int32(m_vcf_hdr, m_bcf_out, "END", &vcf_end_pos, 1);
  }
  //Compute median for all INFO fields
  for(auto i=0u;i<m_INFO_fields_vec.size();++i)
  {
    auto& curr_tuple = m_INFO_fields_vec[i];
    auto known_field_enum = BCF_INFO_GET_KNOWN_FIELD_ENUM(curr_tuple);
    assert(known_field_enum < g_known_variant_field_names.size());
    auto variant_type_enum = BCF_INFO_GET_VARIANT_FIELD_TYPE_ENUM(curr_tuple);
    //valid field handler
    assert(variant_type_enum < m_field_handlers.size() && m_field_handlers[variant_type_enum].get());
    //Just need a 4-byte value, the contents could be a float or int (determined by the templated median function)
    int32_t median;
    auto valid_median_found = m_field_handlers[variant_type_enum]->get_valid_median(copy_variant, *m_query_config,
        m_query_config->get_query_idx_for_known_field_enum(known_field_enum), reinterpret_cast<void*>(&median));
    if(valid_median_found)
      bcf_update_info(m_vcf_hdr, m_bcf_out, g_known_variant_field_names[known_field_enum].c_str(), &median, 1, BCF_INFO_GET_BCF_HT_TYPE(curr_tuple));
  }
}

void BroadCombinedGVCFOperator::handle_FORMAT_fields()
{
  auto& copy_variant = m_variants[0];
  //For weird DP field handling
  auto valid_DP_found = false;
  auto valid_MIN_DP_found = false; 
  m_MIN_DP_vector.resize(copy_variant.get_num_calls()); //will do nothing after the first resize
  auto valid_DP_FORMAT_found = false; 
  m_DP_FORMAT_vector.resize(copy_variant.get_num_calls());
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
    auto variant_type_enum = BCF_FORMAT_GET_VARIANT_FIELD_TYPE_ENUM(curr_tuple);
    //valid field handler
    assert(variant_type_enum < m_field_handlers.size() && m_field_handlers[variant_type_enum].get()); 
    auto valid_field_found = m_field_handlers[variant_type_enum]->collect_and_extend_fields(copy_variant, *m_query_config,
        m_query_config->get_query_idx_for_known_field_enum(known_field_enum), &ptr, num_elements);
    if(valid_field_found)
    {
      auto j=0u;
      auto do_insert = true;    //by default, insert into VCF record
      switch(known_field_enum)
      {
        case GVCF_GT_IDX: //GT field is a pita
          int_vec = const_cast<int*>(reinterpret_cast<const int*>(ptr));
          for(j=0u;j<num_elements;++j)
            int_vec[j] = (int_vec[j] == get_bcf_missing_value<int>()) ? bcf_gt_missing : bcf_gt_unphased(int_vec[j]);
          break;
        case GVCF_GQ_IDX:
          do_insert = m_should_add_GQ_field;
          break;
        case GVCF_MIN_DP_IDX: //simply copy over min-dp values
          memcpy(&(m_MIN_DP_vector[0]), ptr, copy_variant.get_num_calls()*sizeof(int));
          valid_MIN_DP_found = true;
          break;
        case GVCF_DP_FORMAT_IDX:
          memcpy(&(m_DP_FORMAT_vector[0]), ptr, copy_variant.get_num_calls()*sizeof(int));
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
    for(auto j=0ull;j<copy_variant.get_num_calls();++j)
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
      sum_INFO_DP += (is_bcf_valid_value<int>(dp_info_val) ? dp_info_val : 0);
    }
    bcf_update_format_int32(m_vcf_hdr, m_bcf_out, "DP", &(m_DP_FORMAT_vector[0]), m_DP_FORMAT_vector.size()); //add DP FORMAT field
    //If valid DP INFO field found, add to INFO list
    if(valid_DP_found)
      bcf_update_info_int32(m_vcf_hdr, m_bcf_out, "DP", &sum_INFO_DP, 1);
  }
}

void BroadCombinedGVCFOperator::operate(Variant& variant, const VariantQueryConfig& query_config)
{
  GA4GHOperator::operate(variant, query_config);
  assert(m_variants.size() == 1u);      //one variant added by GA4GHOperator::operate, re-orders PL, AD etc
  auto& copy_variant = m_variants[0];
  //Moved to new contig
  if(copy_variant.get_column_begin() >= m_next_contig_begin_position)
    switch_contig();
  //clear out
  bcf_clear(m_bcf_out);
  //position
  m_bcf_out->rid = m_curr_contig_hdr_idx;
  m_bcf_out->pos = copy_variant.get_column_begin() - m_curr_contig_begin_position;
  //Update alleles
  auto& ref_allele = dynamic_cast<VariantFieldString*>(copy_variant.get_common_field(0u).get())->get();
  if(ref_allele.length() == 1u && ref_allele[0] == 'N')
    ref_allele[0] = m_vcf_adapter->get_reference_base_at_position(m_curr_contig_name.c_str(), m_bcf_out->pos);
  const auto& alt_alleles = dynamic_cast<VariantFieldALTData*>(copy_variant.get_common_field(1u).get())->get();
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
  m_should_add_GQ_field = (m_NON_REF_exists && alt_alleles.size() == 1u);
  //INFO fields
  handle_INFO_fields();
  //FORMAT fields
  handle_FORMAT_fields();
  m_vcf_adapter->print_bcf_line(m_bcf_out);
  //Last line in this function always
  m_variants.clear();
}

void BroadCombinedGVCFOperator::switch_contig()
{
  m_curr_contig_name = std::move(m_next_contig_name);
  m_curr_contig_begin_position = m_next_contig_begin_position;
  m_curr_contig_hdr_idx = bcf_hdr_id2int(m_vcf_hdr, BCF_DT_CTG, m_curr_contig_name.c_str());
  auto next_contig_flag = m_vcf_adapter->get_next_contig_location(m_next_contig_begin_position, m_next_contig_name, m_next_contig_begin_position);
  assert(next_contig_flag);
}

#endif //ifdef HTSDIR
