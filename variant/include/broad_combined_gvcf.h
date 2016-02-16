#ifndef BROAD_COMBINED_GVCF_OPERATOR_H
#define BROAD_COMBINED_GVCF_OPERATOR_H

#ifdef HTSDIR

#include "variant_operations.h"
#include "vcf_adapter.h"
#include "vid_mapper.h"

typedef std::tuple<unsigned, unsigned, unsigned> INFO_tuple_type;
typedef std::tuple<unsigned, unsigned, unsigned> FORMAT_tuple_type;

//Exceptions thrown 
class BroadCombinedGVCFException : public std::exception {
  public:
    BroadCombinedGVCFException(const std::string m="") : msg_("Broad combine GVCFs exception : "+m) { ; }
    ~BroadCombinedGVCFException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};


/*
 * Operator to produce the combined GVCF that Broad expects
 */
class BroadCombinedGVCFOperator : public GA4GHOperator
{
  public:
    BroadCombinedGVCFOperator(VCFAdapter& vcf_adapter, const VidMapper& id_mapper, const VariantQueryConfig& query_config);
    ~BroadCombinedGVCFOperator()
    {
      bcf_destroy(m_bcf_out);
      clear();
    }
    void clear();
    void switch_contig();
    virtual void operate(Variant& variant, const VariantQueryConfig& query_config);
    void handle_INFO_fields(const Variant& variant);
    void handle_FORMAT_fields(const Variant& variant);
    void handle_deletions(Variant& variant, const VariantQueryConfig& query_config);
  private:
    const VariantQueryConfig* m_query_config;
    VCFAdapter* m_vcf_adapter;
    const VidMapper* m_vid_mapper;
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
    bool m_should_add_GQ_field;
    //INFO fields enum vector
    std::vector<INFO_tuple_type> m_INFO_fields_vec;
    std::vector<FORMAT_tuple_type> m_FORMAT_fields_vec;
    //MIN_DP values
    std::vector<int> m_MIN_DP_vector;
    //DP_FORMAT values
    std::vector<int> m_DP_FORMAT_vector;
    //Used for handling deletions - remapping PL/AD where a deletion is replaced with *
    CombineAllelesLUT m_reduced_alleles_LUT;
    //vector of field pointers used for handling remapped fields when dealing with spanning deletions
    //avoids re-allocation overhead
    std::vector<std::unique_ptr<VariantFieldBase>> m_spanning_deletions_remapped_fields;
    //Allowed bases
    static const std::unordered_set<char> m_legal_bases;
};

#endif //ifdef HTSDIR

#endif
