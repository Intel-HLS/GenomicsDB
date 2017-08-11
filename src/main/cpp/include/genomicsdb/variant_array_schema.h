/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

#ifndef VARIANT_ARRAY_SCHEMA_H
#define VARIANT_ARRAY_SCHEMA_H

#include "headers.h"
#include "tiledb.h"
#include "gt_common.h"

//Exceptions thrown 
class VariantArraySchemaException : public std::exception {
  public:
    VariantArraySchemaException(const std::string m="") : msg_("VariantArraySchemaException exception : "+m) { ; }
    ~VariantArraySchemaException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

enum VariantArraySchemaFixedFieldsEnum
{
  VARIANT_ARRAY_SCHEMA_END_IDX=0u,
  VARIANT_ARRAY_SCHEMA_REF_IDX,
  VARIANT_ARRAY_SCHEMA_ALT_IDX,
  VARIANT_ARRAY_SCHEMA_QUAL_IDX,
  VARIANT_ARRAY_SCHEMA_FILTER_IDX
};

class AttributeInfo
{
  public:
    AttributeInfo()
      : m_type(typeid(void))
    {
      m_idx = -1;
      m_length = -1;
      m_compression_type = -1;
      m_element_size = UNDEFINED_UINT64_T_VALUE;
    }
    int m_idx;
    int m_length;
    int m_compression_type;
    std::string m_name;
    std::type_index m_type;
    size_t m_element_size;
};

class VariantArraySchema
{
  public:
    /*
     * Empty variant schema
     */
    VariantArraySchema()
      : m_dim_type(std::type_index(typeid(int64_t))) { }
    /*
     * Wrapper around ArraySchema for irregular tiles
     */ 
    VariantArraySchema(const std::string& array_name,
        const std::vector<std::string>& attribute_names,
        const std::vector<std::string>& dim_names,
        const std::vector<std::pair<int64_t, int64_t> >& dim_domains,
        const std::vector<std::type_index>& types,
        const std::vector<int>& val_num, 
        const std::vector<int> compression,
        int cell_order = TILEDB_COL_MAJOR);
    /*
     * Access functions
     * */
    inline const std::string& array_name() const { return m_array_name; }
    inline size_t attribute_num() const { return m_attributes_vector.size(); }
    inline const std::string& attribute_name(int idx) const
    {
      assert(static_cast<size_t>(idx) < m_attributes_vector.size());
      return m_attributes_vector[idx].m_name;
    }
    inline const std::type_index& type(int idx) const
    {
      assert(static_cast<size_t>(idx) < m_attributes_vector.size());
      return m_attributes_vector[idx].m_type;
    }
    inline const int val_num(int idx) const
    {
      assert(static_cast<size_t>(idx) < m_attributes_vector.size());
      return m_attributes_vector[idx].m_length;
    }
    inline const bool is_variable_length_field(const int idx) const
    {
      return (val_num(idx) == TILEDB_VAR_NUM);
    }
    inline const size_t element_size(const int idx) const
    {
      assert(static_cast<size_t>(idx) < m_attributes_vector.size());
      assert(m_attributes_vector[idx].m_element_size != UNDEFINED_UINT64_T_VALUE);
      return m_attributes_vector[idx].m_element_size;
    }
    inline const int compression(int idx) const
    {
      assert(static_cast<size_t>(idx) < m_attributes_vector.size());
      return m_attributes_vector[idx].m_compression_type;
    }
    inline const std::vector<std::pair<int64_t,int64_t>>& dim_domains() const { return m_dim_domains; }
    inline const std::vector<std::string>& dim_names() const { return m_dim_names; }
    inline const std::type_index& dim_type() const { return m_dim_type; }
    inline const int dim_compression_type() const { return m_dim_compression_type; }
    inline const size_t dim_size_in_bytes() const { return m_dim_size_in_bytes; }
  private:
    std::string m_array_name;
    int m_cell_order;
    //Attributes info
    std::vector<AttributeInfo> m_attributes_vector;
    std::unordered_map<std::string, size_t> m_attribute_name_to_idx;
    //Dimensions info
    std::vector<std::pair<int64_t, int64_t>> m_dim_domains;
    std::vector<std::string> m_dim_names;
    std::type_index m_dim_type;
    int m_dim_compression_type;
    size_t m_dim_size_in_bytes;
};

#endif
