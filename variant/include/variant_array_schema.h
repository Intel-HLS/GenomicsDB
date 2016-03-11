#ifndef VARIANT_ARRAY_SCHEMA_H
#define VARIANT_ARRAY_SCHEMA_H

#include "headers.h"
#include "array_schema.h"

class VariantArraySchema
{
  public:
    /*
     * Wrapper around ArraySchema for irregular tiles
     */ 
    VariantArraySchema(const std::string& array_name,
        const std::vector<std::string>& attribute_names,
        const std::vector<std::string>& dim_names,
        const std::vector<std::pair<double, double> >& dim_domains,
        const std::vector<const std::type_info*>& types,
        const std::vector<int>& val_num, 
        const std::vector<CompressionType> compression,
        ArraySchema::CellOrder cell_order = ArraySchema::CO_COLUMN_MAJOR,
        int64_t capacity = AS_CAPACITY,
        int consolidation_step = AS_CONSOLIDATION_STEP)
    {
      m_array_schema = new ArraySchema(array_name, attribute_names, dim_names, dim_domains,
          types, val_num, compression, cell_order, capacity, consolidation_step);
      m_owns_schema = true;
    }
    /*
     * Un-owned array schema constructor
     */
    VariantArraySchema(const ArraySchema* array_schema)
    {
      m_array_schema = array_schema;
      m_owns_schema = false;
    }
    /*
     * "Copy" constructor
     */
    VariantArraySchema(const VariantArraySchema& other)
    {
      m_array_schema = other.m_array_schema;
      m_owns_schema = false;
    }
    /*
     * "Assignment" operator
     */
    VariantArraySchema& operator=(const VariantArraySchema& other)
    {
      if(m_owns_schema && m_array_schema)
        delete m_array_schema;
      m_array_schema = other.m_array_schema;
      m_owns_schema = false;
      return *this;
    }
    /*
     * Destructor
     */
    ~VariantArraySchema()
    {
      if(m_owns_schema && m_array_schema)
        delete m_array_schema;
      m_array_schema = 0;
      m_owns_schema = false;
    }
    inline const ArraySchema* get_array_schema() const { return m_array_schema; }
    inline int attribute_num() const { return m_array_schema->attribute_num(); }
    inline const std::string& attribute_name(int idx) const { return m_array_schema->attribute_name(idx); }
    inline const std::type_info* type(int idx) const { return m_array_schema->type(idx); }
    inline const int val_num(int idx) const { return m_array_schema->val_num(idx); }
    inline const std::vector<std::pair<double,double>>& dim_domains() const { return m_array_schema->dim_domains(); }
    inline const std::string& array_name() const { return m_array_schema->array_name(); }
  private:
    const ArraySchema* m_array_schema;
    bool m_owns_schema;
};

#endif
