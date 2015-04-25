#ifndef GT_COMMON_H
#define GT_COMMON_H

#include <stdio.h>
#include <typeinfo>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include <assert.h>
#include <math.h>
#include <queue>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <string>
#include "profiling.h"

#ifndef HTSDIR
#include "vcf.h"
#else
#include "htslib/vcf.h"
#endif

#define CHECK_MISSING_SAMPLE_GIVEN_REF(REF) (((REF).size() == 0) || ((REF)[0] == '$'))
#define CHECK_UNINITIALIZED_SAMPLE_GIVEN_REF(REF) ((REF).size() == 0 || ((REF) == ""))
#define CHECK_IN_THE_MIDDLE_REF(REF) ((REF)[0] == 'N')

inline bool IS_NON_REF_ALLELE(const std::string& allele)
{
  return allele.length() > 0 && ((allele)[0] == '&');
}

inline bool IS_NON_REF_ALLELE(const char allele_char)
{
  return allele_char == '&';
}

#define TILEDB_NON_REF_VARIANT_REPRESENTATION "&"

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
  GVCF_SB_1_IDX,
  GVCF_SB_2_IDX,
  GVCF_SB_3_IDX,
  GVCF_SB_4_IDX,
  GVCF_AD_IDX,
  GVCF_PL_IDX,
  GVCF_AF_IDX,
  GVCF_AN_IDX,
  GVCF_AC_IDX,
  GVCF_NULL_IDX,
  GVCF_OFFSETS_IDX,
  GVCF_COORDINATES_IDX,
  GVCF_NUM_KNOWN_FIELDS
};

#endif
