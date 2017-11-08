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

#ifndef KNOWN_FIELD_INFO_H
#define KNOWN_FIELD_INFO_H

#include "gt_common.h"
#include <memory>

//Should be identical to the vector m_known_variant_field_names, see file known_field_info.cc
enum KnownVariantFieldsEnum
{
  GVCF_END_IDX = 0,
  GVCF_REF_IDX,
  GVCF_ALT_IDX,
  GVCF_QUAL_IDX,
  GVCF_FILTER_IDX,
  GVCF_BASEQRANKSUM_IDX,
  GVCF_CLIPPINGRANKSUM_IDX,
  GVCF_MQRANKSUM_IDX,
  GVCF_READPOSRANKSUM_IDX,
  GVCF_DP_IDX,
  GVCF_MQ_IDX,
  GVCF_RAW_MQ_IDX,
  GVCF_MQ0_IDX,
  GVCF_DP_FORMAT_IDX,
  GVCF_MIN_DP_IDX,
  GVCF_GQ_IDX,
  GVCF_SB_IDX,
  GVCF_AD_IDX,
  GVCF_PL_IDX,
  GVCF_AF_IDX,
  GVCF_AN_IDX,
  GVCF_AC_IDX,
  GVCF_GT_IDX,
  GVCF_PS_IDX,
  GVCF_PGT_IDX,
  GVCF_PID_IDX,
  GVCF_EXCESS_HET,
  GVCF_ID_IDX,
  GVCF_NUM_KNOWN_FIELDS
};

//All known field names specific to variant data
extern std::vector<std::string> g_known_variant_field_names;
//Mapping from field name to enum idx
extern std::unordered_map<std::string, unsigned> g_known_variant_field_name_to_enum;

//Known fields exception
class KnownFieldInfoException : public std::exception {
  public:
    KnownFieldInfoException(const std::string m="") : msg_("KnownFieldInfoException : "+m) { ; }
    ~KnownFieldInfoException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class VariantFieldCreatorBase;
/*
 * Class that stores info about some of the known fields
 */
class KnownFieldInfo
{
  friend class KnownFieldInitializer;
  public:
    KnownFieldInfo();
  private:
    unsigned m_length_descriptor;
    unsigned m_num_elements;
    std::shared_ptr<VariantFieldCreatorBase> m_field_creator;
    int m_VCF_field_combine_operation;
  public:
    inline bool is_length_allele_dependent() const
    {
      unsigned length_descriptor = m_length_descriptor;
      return (length_descriptor == BCF_VL_A || length_descriptor == BCF_VL_R || length_descriptor == BCF_VL_G);
    }
    inline unsigned get_length_descriptor() const { return m_length_descriptor; }
    inline bool is_length_genotype_dependent() const { return m_length_descriptor == BCF_VL_G; }
    inline bool is_length_only_ALT_alleles_dependent() const { return m_length_descriptor == BCF_VL_A; }
    unsigned get_num_elements_for_known_field_enum(unsigned num_ALT_alleles, unsigned ploidy) const;
    inline int get_VCF_field_combine_operation() const { return m_VCF_field_combine_operation; }
    /*
     * Static functions that access the global vector specified below to get info
     */
    inline static bool is_length_descriptor_genotype_dependent(unsigned length_descriptor) { return length_descriptor == BCF_VL_G; }
    inline static bool is_length_descriptor_only_ALT_alleles_dependent(unsigned length_descriptor) { return length_descriptor == BCF_VL_A; }
    inline static bool is_length_descriptor_all_alleles_dependent(unsigned length_descriptor) { return length_descriptor == BCF_VL_R; }
    inline static bool is_length_descriptor_ploidy_dependent(unsigned length_descriptor)
    { return ((length_descriptor == BCF_VL_P) || (length_descriptor == BCF_VL_Phased_Ploidy)); }
    /*
     * Given a field name, checks for m_known_variant_field_name_to_enum to see if valid entry exists.
     * If yes, fills known_field_enum and returns true
     * Else returns false. Leaves known_field_enum untouched.
     */
    static bool get_known_field_enum_for_name(const std::string& field_name, unsigned& known_field_enum);
    /*
     * Get name for known field enum
     */
    static std::string get_known_field_name_for_enum(unsigned known_field_enum);
    /*
     * Check whether the known field requires a special creator
     */
    static bool requires_special_creator(unsigned enumIdx);
    static const std::shared_ptr<VariantFieldCreatorBase>& get_field_creator(unsigned enumIdx);
    /*
     * Function that determines whether length of the field is dependent on the #alleles 
     */
    static bool is_length_allele_dependent(unsigned enumIdx);
    /*
     * Function that determines whether length descriptor is dependent on the #alleles 
     */
    static bool is_length_descriptor_allele_dependent(unsigned length_descriptor)
    {
      return (length_descriptor == BCF_VL_A || length_descriptor == BCF_VL_R || length_descriptor == BCF_VL_G);
    }
    /*
     * Return ploidy given length descriptor and the number of elements in the cell for the GT field
     */
    static inline unsigned get_ploidy(const unsigned length_descriptor, const unsigned num_elements)
    {
      switch(length_descriptor)
      {
        case BCF_VL_P:
          return num_elements;
        case BCF_VL_Phased_Ploidy:
          return ((num_elements+1u) >> 1u); //Eg. 0/1 becomes [0,0,1], 0|2/1 becomes [0,1,2,0,1]
        default:
          throw KnownFieldInfoException(std::string("Unknown length descriptor for GT field ")
              + std::to_string(length_descriptor));
          return 0u;
      }
    }
    /*
     * Given a length descriptor, get #elements
     */
    static unsigned get_num_elements_given_length_descriptor(unsigned length_descriptor,
        unsigned num_ALT_alleles, unsigned ploidy, unsigned num_elements);
    /*
     * Get #possible genotypes given #ALT alleles and ploidy
     */
    static unsigned get_number_of_genotypes(const unsigned num_ALT_alleles, const unsigned ploidy);
    /*
     * Function that determines whether length of the field is dependent on the #genotypes
     */
    static bool is_length_genotype_dependent(unsigned enumIdx);
    /*
     * Function that determines whether length of the field is dependent only on the #alt alleles
     */
    static bool is_length_only_ALT_alleles_dependent(unsigned enumIdx);
    /*
     * Functions that determine number of elements for known fields
     */

    static unsigned get_num_elements_for_known_field_enum(unsigned known_field_enum,
        unsigned num_ALT_alleles, unsigned ploidy);
    /*
     * Returns BCF length descriptor, BCF_VL_*
     */
    static unsigned get_length_descriptor_for_known_field_enum(unsigned known_field_enum);
    /*
     * INFO field combine operation
     */ 
    static int get_VCF_field_combine_operation_for_known_field_enum(unsigned known_field_enum);
};
/*
 * Vector that stores information about the known fields - length, Factory methods etc
 */
extern std::vector<KnownFieldInfo> g_known_field_enum_to_info;
/*
 * Class whose sole purpose is to initialize KnownFieldInfo etc
 * Only a single instance of this class should exist in the whole program - the global variable listed below
 * The constructor initializes all KnownFieldInfo
 */
class KnownFieldInitializer
{
  public:
    KnownFieldInitializer();
  private:
    void initialize_length_descriptor(unsigned idx) const;
    void initialize_INFO_combine_operation(unsigned idx) const;
};
extern KnownFieldInitializer g_known_field_initializer;

class VariantUtils
{
  public:
    static bool contains_deletion(const std::string& REF, const std::vector<std::string>& ALT_vec);
    inline static bool is_deletion(const std::string& REF, const std::string& alt_allele)
    {
      auto REF_length = REF.length();
      return (REF_length > 1u) && !IS_NON_REF_ALLELE(alt_allele) && (alt_allele.length() < REF_length);
    }
    inline static bool is_symbolic_allele(const std::string& allele)
    {
      return IS_NON_REF_ALLELE(allele)
        || (allele == g_vcf_SPANNING_DELETION)
        || (
            allele.length() > 0u &&
            (
             (allele[0] == '<' && allele[allele.length()-1u] == '>') || 
             (allele.find_first_of('[') != std::string::npos || allele.find_first_of(']') != std::string::npos)
            ) 
           );
    }
    inline static bool is_reference_block(const std::string& REF, const std::vector<std::string>& ALT_vec)
    {
      return (REF.length() == 1u && ALT_vec.size() == 1u && IS_NON_REF_ALLELE(ALT_vec[0]));
    }
};

#endif
