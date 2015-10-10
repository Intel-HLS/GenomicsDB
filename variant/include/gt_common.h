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

#endif
