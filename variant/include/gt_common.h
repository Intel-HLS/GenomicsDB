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

enum GTSchemaVersionEnum
{
    GT_SCHEMA_V1=0,
    GT_SCHEMA_V2
};

#ifdef DO_PROFILING
extern uint64_t g_num_disk_loads;
extern uint64_t g_num_cached_loads;
extern uint64_t g_coords_num_disk_loads;
extern uint64_t g_coords_num_cached_loads;
extern uint64_t g_total_num_tiles_loaded;
extern std::vector<uint64_t> g_num_tiles_loaded;
extern std::vector<uint64_t> g_num_segments_loaded;
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
