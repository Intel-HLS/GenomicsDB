#ifndef GT_COMMON_H
#define GT_COMMON_H

#include "headers.h"
#include "profiling.h"

#include "vcf.h"

#define CHECK_MISSING_SAMPLE_GIVEN_REF(REF) (((REF).size() == 0) || ((REF)[0] == '$'))
#define CHECK_UNINITIALIZED_SAMPLE_GIVEN_REF(REF) ((REF).size() == 0 || ((REF) == ""))
#define CHECK_IN_THE_MIDDLE_REF(REF) ((REF)[0] == 'N')
extern std::string g_vcf_NON_REF;
inline bool IS_NON_REF_ALLELE(const std::string& allele)
{
  return allele.length() > 0 && ((allele)[0] == '&');
}

inline bool IS_NON_REF_ALLELE(const char allele_char)
{
  return allele_char == '&';
}

#define TILEDB_NON_REF_VARIANT_REPRESENTATION "&"
#define TILEDB_ALT_ALLELE_SEPARATOR "|"

//Should be identical to the vector m_known_variant_field_names, see file query_variants.cc
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
  GVCF_MQ0_IDX,
  GVCF_DP_FMT_IDX,
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
  GVCF_NUM_KNOWN_FIELDS
};

//All known field names specific to variant data
extern std::vector<std::string> g_known_variant_field_names;
//Mapping from field name to enum idx
extern std::unordered_map<std::string, unsigned> g_known_variant_field_name_to_enum;

#endif
