#ifndef VARIANT_ARRAY_SCHEMA_H
#define VARIANT_ARRAY_SCHEMA_H

#include "headers.h"
#include "c_api.h"

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

class AttributeInfo
{
  public:
    AttributeInfo() : m_type(typeid(void)) {}
    int m_idx;
    int m_length;
    int m_compression_type;
    std::string m_name;
    std::type_index m_type;
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
