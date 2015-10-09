#ifdef HTSDIR

#include "broad_combined_gvcf.h"

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
  m_INFO_fields_vec = std::move(std::vector<unsigned>
      {
      GVCF_BASEQRANKSUM_IDX,
      GVCF_MQ_IDX,
      GVCF_MQ0_IDX,
      GVCF_CLIPPINGRANKSUM_IDX,
      GVCF_MQRANKSUM_IDX,
      GVCF_READPOSRANKSUM_IDX,
      GVCF_DP_IDX
      });
  m_FORMAT_fields_vec = std::move(std::vector<unsigned>
      {
      GVCF_GT_IDX,
      GVCF_DP_FMT_IDX,
      GVCF_GQ_IDX,
      GVCF_MIN_DP_IDX,
      GVCF_AD_IDX,
      GVCF_PL_IDX,
      GVCF_SB_IDX
      });
  //TODO: add samples to header
}

void BroadCombinedGVCFOperator::clear()
{
  m_curr_contig_name.clear();
  m_next_contig_name.clear();
  m_alleles_pointer_buffer.clear();
  m_INFO_fields_vec.clear();
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
  auto total_num_merged_alleles = m_merged_alt_alleles.size() + 1u;      //+1 for REF
  if(total_num_merged_alleles > m_alleles_pointer_buffer.size())
    m_alleles_pointer_buffer.resize(total_num_merged_alleles);
  //REF
  m_alleles_pointer_buffer[0] = m_merged_reference_allele.c_str();
  //ALT
  for(auto i=1u;i<total_num_merged_alleles;++i)
    m_alleles_pointer_buffer[i] = m_merged_alt_alleles[i-1u].c_str();
  bcf_update_alleles(m_vcf_hdr, m_bcf_out, &(m_alleles_pointer_buffer[0]), total_num_merged_alleles);
  handle_INFO_fields();
  m_vcf_adapter->print_bcf_line(m_bcf_out);
  //Last line always
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
