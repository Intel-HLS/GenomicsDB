#include "variant_field_data.h"

std::string g_vcf_NON_REF="<NON_REF>";

fi_union g_tiledb_null_float = { .f = TILEDB_EMPTY_FLOAT32 };
di_union g_tiledb_null_double = { .d = TILEDB_EMPTY_FLOAT64 };

std::unordered_map<std::type_index, VariantFieldTypeEnum> g_variant_field_type_index_to_enum = 
std::unordered_map<std::type_index, VariantFieldTypeEnum>{
  { std::type_index(typeid(void)), VARIANT_FIELD_VOID },
  { std::type_index(typeid(int)), VARIANT_FIELD_INT },
  { std::type_index(typeid(int64_t)), VARIANT_FIELD_INT64_T },
  { std::type_index(typeid(unsigned)), VARIANT_FIELD_UNSIGNED },
  { std::type_index(typeid(uint64_t)), VARIANT_FIELD_UINT64_T },
  { std::type_index(typeid(float)), VARIANT_FIELD_FLOAT },
  { std::type_index(typeid(double)), VARIANT_FIELD_DOUBLE },
  { std::type_index(typeid(std::string)), VARIANT_FIELD_STRING },
  { std::type_index(typeid(char)), VARIANT_FIELD_CHAR }
};

std::unordered_map<std::type_index, int> g_variant_field_type_index_to_tiledb_type = 
std::unordered_map<std::type_index, int>{
  { std::type_index(typeid(int)), TILEDB_INT32 },
  { std::type_index(typeid(int64_t)), TILEDB_INT64 },
  { std::type_index(typeid(unsigned)), TILEDB_INT32 },
  { std::type_index(typeid(uint64_t)), TILEDB_INT64 },
  { std::type_index(typeid(float)), TILEDB_FLOAT32 },
  { std::type_index(typeid(double)), TILEDB_FLOAT64 },
  { std::type_index(typeid(char)), TILEDB_CHAR }
};

std::vector<std::type_index> g_tiledb_type_to_variant_field_type_index =
{
  std::type_index(typeid(int)),
  std::type_index(typeid(int64_t)),
  std::type_index(typeid(float)),
  std::type_index(typeid(double)),
  std::type_index(typeid(char))
};

size_t VariantFieldTypeUtil::size(const VariantFieldTypeEnum type_enum)
{
  switch(type_enum)
  {
    case VARIANT_FIELD_INT:
      return sizeof(int);
    case VARIANT_FIELD_INT64_T:
      return sizeof(int64_t);
    case VARIANT_FIELD_UNSIGNED:
      return sizeof(unsigned);
    case VARIANT_FIELD_UINT64_T:
      return sizeof(uint64_t);
    case VARIANT_FIELD_FLOAT:
      return sizeof(float);
    case VARIANT_FIELD_DOUBLE:
      return sizeof(double);
    case VARIANT_FIELD_STRING:
      return sizeof(std::string);
    case VARIANT_FIELD_CHAR:
      return sizeof(char);
    default:
      return 0;
  }
  return 0;
}
