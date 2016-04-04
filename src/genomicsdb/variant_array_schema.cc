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

#include "variant_array_schema.h"
#include "variant_field_data.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw VariantArraySchemaException(#X);

VariantArraySchema::VariantArraySchema(const std::string& array_name,
        const std::vector<std::string>& attribute_names,
        const std::vector<std::string>& dim_names,
        const std::vector<std::pair<int64_t, int64_t> >& dim_domains,
        const std::vector<std::type_index>& types,
        const std::vector<int>& val_num, 
        const std::vector<int> compression,
        int cell_order)
  : m_dim_type(typeid(int64_t))
{
  m_array_name = array_name;
  m_cell_order = cell_order;
  VERIFY_OR_THROW(attribute_names.size() == val_num.size());
  VERIFY_OR_THROW(attribute_names.size()+1u == types.size() &&
      compression.size() == types.size() &&
      "Last element of types and compression vectors must specify type and compression of co-ordinates");
  VERIFY_OR_THROW(dim_names.size() == dim_domains.size());
  m_attributes_vector.resize(attribute_names.size());
  for(auto i=0u;i<attribute_names.size();++i)
  {
    m_attribute_name_to_idx[attribute_names[i]] = i;
    auto& curr_elem = m_attributes_vector[i];
    curr_elem.m_idx = i;
    curr_elem.m_length = val_num[i];
    curr_elem.m_compression_type = compression[i];
    curr_elem.m_name = attribute_names[i];
    curr_elem.m_type = types[i];
    curr_elem.m_element_size = VariantFieldTypeUtil::size(types[i]);
  }
  //Co-ordinates
  m_dim_names = dim_names;
  m_dim_domains = dim_domains;
  m_dim_type = std::type_index(types[types.size()-1u]);
  m_dim_compression_type = compression[compression.size()-1u];
  m_dim_size_in_bytes = m_dim_names.size()*VariantFieldTypeUtil::size(m_dim_type);
}

