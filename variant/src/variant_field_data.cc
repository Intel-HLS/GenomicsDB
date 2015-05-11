#include "variant_field_data.h"

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
