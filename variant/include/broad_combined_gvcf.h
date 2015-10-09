#ifndef BROAD_COMBINED_GVCF_OPERATOR_H
#define BROAD_COMBINED_GVCF_OPERATOR_H

#ifdef HTSDIR

#include "variant_operations.h"
#include "vcf_adapter.h"

/*
 * Operator to produce the combined GVCF that Broad expects
 */
class BroadCombinedGVCFOperator : public GA4GHOperator
{
  public:
    BroadCombinedGVCFOperator(VCFAdapter& vcf_adapter, const VariantQueryConfig& query_config);
    ~BroadCombinedGVCFOperator()
    {
      bcf_destroy(m_bcf_out);
      clear();
    }
    void clear();
    void switch_contig();
    virtual void operate(Variant& variant, const VariantQueryConfig& query_config);
    void handle_INFO_fields();
  private:
    const VariantQueryConfig* m_query_config;
    VCFAdapter* m_vcf_adapter;
    bcf1_t* m_bcf_out;
    bcf_hdr_t* m_vcf_hdr;
    //Contig info
    std::string m_curr_contig_name;
    int64_t m_curr_contig_begin_position;
    int m_curr_contig_hdr_idx;
    std::string m_next_contig_name;
    int64_t m_next_contig_begin_position;
    //alleles pointers buffer
    std::vector<const char*> m_alleles_pointer_buffer;
    //INFO fields enum vector
    std::vector<unsigned> m_INFO_fields_vec;
    std::vector<unsigned> m_FORMAT_fields_vec;
};

#endif //ifdef HTSDIR

#endif
