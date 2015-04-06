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

enum GTSchemaVersionEnum
{
    GT_SCHEMA_V1=0,
    GT_SCHEMA_V2
};

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
