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
#include <algorithm>
#include <string>

enum GVCFAttributesEnum
{
  GVCF_END_IDX=0,
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
  GVCF_NULL_IDX,
  GVCF_OFFSETS_IDX,
  GVCF_COORDINATES_IDX
};

enum GTAttributeIdx
{
  GT_END_IDX=0,
  GT_REF_IDX,
  GT_ALT_IDX,
  GT_PL_IDX,
  GT_NULL_IDX,
  GT_OFFSETS_IDX,
  GT_COORDINATES_IDX
};

#ifdef DO_PROFILING
extern uint64_t g_num_disk_loads;
extern uint64_t g_num_cached_loads;
extern uint64_t g_coords_num_disk_loads;
extern uint64_t g_coords_num_cached_loads;
extern uint64_t g_total_num_tiles_loaded;
extern uint64_t g_num_tiles_loaded[GVCF_COORDINATES_IDX+1];
extern uint64_t g_num_segments_loaded[GVCF_COORDINATES_IDX+1];
#endif

#ifndef HTSDIR
#include "vcf.h"
#else
#include "htslib/vcf.h"
#endif

#define CHECK_MISSING_SAMPLE_GIVEN_REF(REF) (((REF).size() == 0) || ((REF)[0] == '$'))
#define CHECK_UNINITIALIZED_SAMPLE_GIVEN_REF(REF) ((REF).size() == 0 || ((REF) == ""))
#define CHECK_IN_THE_MIDDLE_REF(REF) ((REF)[0] == 'N')
#define IS_NON_REF_ALLELE(allele) ((allele)[0] == '&')


#endif
