#ifndef GT_COMMON_H
#define GT_COMMON_H

#include "headers.h"
#include "profiling.h"

#include "vcf.h"

#define CHECK_MISSING_SAMPLE_GIVEN_REF(REF) (((REF).size() == 0) || ((REF)[0] == '$'))
#define CHECK_UNINITIALIZED_SAMPLE_GIVEN_REF(REF) ((REF).size() == 0 || ((REF) == ""))
#define CHECK_IN_THE_MIDDLE_REF(REF) ((REF)[0] == 'N')
extern std::string g_vcf_NON_REF;
extern std::string g_vcf_SPANNING_DELETION;
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

#define UNDEFINED_ATTRIBUTE_IDX_VALUE 0xFFFFFFFFu
#define UNDEFINED_NUM_ROWS_VALUE 0xFFFFFFFFFFFFFFFFull

enum VariantFieldTypeEnum
{
  VARIANT_FIELD_VOID=0,
  VARIANT_FIELD_INT,
  VARIANT_FIELD_INT64_T,
  VARIANT_FIELD_UNSIGNED,
  VARIANT_FIELD_UINT64_T,
  VARIANT_FIELD_FLOAT,
  VARIANT_FIELD_DOUBLE,
  VARIANT_FIELD_STRING,
  VARIANT_FIELD_CHAR,
  VARIANT_FIELD_NUM_TYPES
};

//Map from type_index to VariantFieldTypeEnum
extern std::unordered_map<std::type_index, VariantFieldTypeEnum> g_variant_field_type_index_to_enum;
extern std::unordered_map<std::type_index, int> g_variant_field_type_index_to_tiledb_type;
extern std::vector<std::type_index> g_tiledb_type_to_variant_field_type_index;

#endif
