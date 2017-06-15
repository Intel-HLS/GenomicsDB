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

#ifndef BROAD_COMBINED_GVCF_OPERATOR_H
#define BROAD_COMBINED_GVCF_OPERATOR_H

#ifdef HTSDIR

#include "variant_operations.h"
#include "vcf_adapter.h"
#include "vid_mapper.h"
#include "timer.h"

//known_field_enum, query_idx, VariantFieldTypeEnum, bcf_ht_type, vcf field name, INFO_field_combine_operation
typedef std::tuple<unsigned, unsigned, VariantFieldTypeEnum, unsigned, std::string, int> INFO_tuple_type;
//known_field_enum, query_idx, VariantFieldTypeEnum, bcf_ht_type, vcf field name
typedef std::tuple<unsigned, unsigned, VariantFieldTypeEnum, unsigned, std::string> FORMAT_tuple_type;

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
    BroadCombinedGVCFOperator(VCFAdapter& vcf_adapter, const VidMapper& id_mapper, const VariantQueryConfig& query_config,
        const unsigned max_diploid_alt_alleles_that_can_be_genotyped=MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED,
        const bool use_missing_values_only_not_vector_end=false);
    virtual ~BroadCombinedGVCFOperator()
    {
      bcf_destroy(m_bcf_out);
      clear();
#ifdef DO_PROFILING
      m_bcf_t_creation_timer.print("bcf_t creation time", std::cerr);
#endif
    }
    void clear();
    void switch_contig();
    virtual void operate(Variant& variant, const VariantQueryConfig& query_config);
    inline bool overflow() const { return m_vcf_adapter->overflow(); }
    bool handle_VCF_field_combine_operation(const Variant& variant,
        const INFO_tuple_type& curr_tuple, void*& result_ptr, unsigned& num_result_elements);
    void handle_INFO_fields(const Variant& variant);
    void handle_FORMAT_fields(const Variant& variant);
    void handle_deletions(Variant& variant, const VariantQueryConfig& query_config);
    void merge_ID_field(const Variant& variant, const unsigned query_idx);
  private:
    bool m_use_missing_values_not_vector_end;
    const VariantQueryConfig* m_query_config;
    VCFAdapter* m_vcf_adapter;
    const VidMapper* m_vid_mapper;
    bcf1_t* m_bcf_out;
    size_t m_bcf_record_size;
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
    //If QUAL combine operation is specified
    INFO_tuple_type m_vcf_qual_tuple;
    //ID VCF field - to avoid repeated dynamic reallocations
    std::string m_ID_value;
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
    std::vector<int> m_spanning_deletion_remapped_GT;
    //Allowed bases
    static const std::unordered_set<char> m_legal_bases;
    //For profiling
    Timer m_bcf_t_creation_timer;
};

#endif //ifdef HTSDIR

#endif
