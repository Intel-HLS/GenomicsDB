/**
 * The MIT License (MIT)
 * Copyright (c) 2016 Intel Corporation
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

#include "genomicsdb_columnar_field.h"
#include "variant_field_data.h"

GenomicsDBColumnarField::GenomicsDBColumnarField(const std::type_index element_type, const int length_descriptor,
    const unsigned fixed_length_field_length)
: m_element_type(element_type), m_length_descriptor(length_descriptor), m_fixed_length_field_num_elements(fixed_length_field_length)
{
  m_element_size = VariantFieldTypeUtil::size(m_element_type);
  m_fixed_length_field_size = m_fixed_length_field_num_elements*m_element_size;
  assign_function_pointers();
}

void GenomicsDBColumnarField::resize(const uint64_t num_cells)
{
  m_valid.resize(num_cells, false);
  if(m_length_descriptor == BCF_VL_FIXED)
    m_data.resize(num_cells*m_element_size*m_fixed_length_field_num_elements);
  else
  {
    m_num_bytes.resize(num_cells, 0ull);
    m_begin_offsets.resize(num_cells+1, 0ull);
  }
}

void GenomicsDBColumnarField::assign_function_pointers()
{
  switch(VariantFieldTypeUtil::get_variant_field_type_enum_for_variant_field_type(m_element_type))
  {
    case VariantFieldTypeEnum::VARIANT_FIELD_INT:
    case VariantFieldTypeEnum::VARIANT_FIELD_UNSIGNED:
    case VariantFieldTypeEnum::VARIANT_FIELD_INT64_T:
    case VariantFieldTypeEnum::VARIANT_FIELD_UINT64_T:
      m_check_validity = GenomicsDBColumnarField::check_validity<int>;
      break;
    case VariantFieldTypeEnum::VARIANT_FIELD_FLOAT:
      m_check_validity = GenomicsDBColumnarField::check_validity<float>;
      break;
    case VariantFieldTypeEnum::VARIANT_FIELD_CHAR:
    case VariantFieldTypeEnum::VARIANT_FIELD_STRING:
      m_check_validity = GenomicsDBColumnarField::check_validity<char>;
      break;
    default:
      throw GenomicsDBColumnarFieldException(std::string("Unhandled type ")+m_element_type.name());
  }
}
