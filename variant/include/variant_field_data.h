#ifndef VARIANT_FIELD_DATA_H
#define VARIANT_FIELD_DATA_H

#include "query_field_data.h"
#include <memory>

/*
 * Base class for variant field data - not sure whether I will add any functionality here
 */
class VariantFieldBase : public QueryFieldData
{
  public:
    VariantFieldBase() 
      : QueryFieldData()
    { ; }
    virtual ~VariantFieldBase() = default;
    virtual void copy_data_from_tile(const Tile&  base_tile, uint64_t element_idx, uint64_t num_elements) = 0;
    virtual void clear() { ; }
    virtual void print(std::ostream& fptr) const { ; }
    /* Get pointer(s) to data with number of elements */
    virtual std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated) = 0;
    /* Create copy and return pointer - avoid using as much as possible*/
    virtual VariantFieldBase* create_copy() const = 0;
};
/*
 * Class that holds single element data 
 */
template<class DataType, class AttributeTileTy=const AttributeTile<DataType>>
class VariantFieldData : public VariantFieldBase
{
  public:
    VariantFieldData()
      : VariantFieldBase()
    { ; }
    virtual ~VariantFieldData() = default;
    virtual void copy_data_from_tile(const Tile& base_tile, uint64_t element_idx, uint64_t num_elements=0)
    {
      const AttributeTileTy& tile = static_cast<const AttributeTileTy&>(base_tile);
      m_data = tile.cell(element_idx);
      /*tile.copy_payload(&m_data, element_idx, 1u);*/
    }
    virtual DataType& get() { return m_data; }
    virtual const DataType& get() const { return m_data; }
    virtual std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated)
    {
      size = 1u;
      *(reinterpret_cast<DataType**>(ptr)) = &m_data;
      allocated = false;
      return std::type_index(typeid(DataType));
    }
    virtual VariantFieldBase* create_copy() const { return new VariantFieldData<DataType, AttributeTileTy>(*this); }
  private:
    DataType m_data;
};
/*
 * Specialize class for string, AttributeTile type is <char> 
 */
template<>
class VariantFieldData<std::string, const AttributeTile<char>> : public VariantFieldBase
{
  public:
    VariantFieldData()
      : VariantFieldBase()
    { ; }
    virtual ~VariantFieldData() = default;
    virtual void clear() { m_data.clear(); }
    virtual void copy_data_from_tile(const Tile& base_tile, uint64_t element_idx, uint64_t num_elements=0)
    {
      const AttributeTile<char>& tile = static_cast<const AttributeTile<char>&>(base_tile);
      /** Assuming the string in the payload is null terminated, else dead.
       * This still copies the data into m_data, effectively strcpy, which hopefully is faster than any
       * custom code we write **/
      const auto* payload_ptr = tile.get_payload_ptr(element_idx);
      assert(element_idx < tile.cell_num());
      auto string_length = strnlen(payload_ptr, tile.cell_num()-element_idx); //do not run off beyond end of tile
      m_data = std::move(std::string(tile.get_payload_ptr(element_idx), string_length));
    }
    virtual std::string& get()  { return m_data; }
    virtual const std::string& get() const { return m_data; }
    virtual void print(std::ostream& fptr) const { fptr << m_data; }
    virtual std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated)
    {
      size = 1u;
      char** data = new char*;
      data[0] = const_cast<char*>(m_data.c_str());
      *(reinterpret_cast<char***>(ptr)) = data;
      allocated = true;
      return std::type_index(typeid(char));
    }
    virtual VariantFieldBase* create_copy() const { return new VariantFieldData<std::string, const AttributeTile<char>>(*this); }
  private:
    std::string m_data;
};
//Assigned name for string type
typedef VariantFieldData<std::string, const AttributeTile<char>> VariantFieldString;
/*
 * Sub-class that holds vector data of basic types - int,float etc
 */
template<class DataType, class AttributeTileTy=const AttributeTile<DataType>>
class VariantFieldPrimitiveVectorData : public VariantFieldBase
{
  public:
    VariantFieldPrimitiveVectorData()
      : VariantFieldBase()
    { clear(); }
    virtual ~VariantFieldPrimitiveVectorData() = default;
    virtual void clear() { m_data.clear(); }
    virtual void copy_data_from_tile(const Tile& base_tile, uint64_t element_idx, uint64_t num_elements)
    {
      const AttributeTileTy& tile = static_cast<const AttributeTileTy&>(base_tile);
      m_data.resize(num_elements);
      tile.copy_payload(&(m_data[0]), element_idx, num_elements);
    }
    virtual std::vector<DataType>& get()  { return m_data; }
    virtual const std::vector<DataType>& get() const { return m_data; }
    virtual void print(std::ostream& fptr) const
    {
      fptr << "[ ";
      for(auto val : m_data)
        fptr << val << ",";
      fptr << "]";
    }
    virtual std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated)
    {
      size = m_data.size();
      *(reinterpret_cast<DataType**>(ptr)) = (size > 0) ? &(m_data[0]) : nullptr;
      allocated = false;
      return std::type_index(typeid(DataType));
    }
    virtual VariantFieldBase* create_copy() const { return new VariantFieldPrimitiveVectorData<DataType, AttributeTileTy>(*this); }
  private:
    std::vector<DataType> m_data;
};
/*
 * Special class for ALT field. ALT field is parsed in a weird way, hence treated in a special manner
 */
class VariantFieldALTData : public VariantFieldBase
{
  public:
    VariantFieldALTData()
      : VariantFieldBase()
    { ; }
    virtual ~VariantFieldALTData() = default;
    virtual void clear()
    {
      for(auto& s : m_data)
        s.clear();
      m_data.clear();
    }
    virtual void copy_data_from_tile(const Tile& base_tile, uint64_t element_idx, uint64_t num_elements=0)
    {
      auto insert_idx = 0u;
      const AttributeTile<char>& tile = static_cast<const AttributeTile<char>&>(base_tile);
      assert(element_idx < tile.cell_num());
      uint64_t offset_idx = element_idx;
      /** Exit if empty string or last element is NON_REF **/
      while(1)
      {
        /* Assuming the string in the payload is null terminated.
         * This still copies the data into m_data, effectively calling strcpy which, hopefully, is faster than any
         * custom code we write **/
        const auto* payload_ptr = tile.get_payload_ptr(offset_idx);
        if(offset_idx >= tile.cell_num())     //end of tile
          break;
        auto string_length = strnlen(payload_ptr, tile.cell_num()-offset_idx); //do not run off beyond end of tile
        std::string last_string = std::move(std::string(payload_ptr, string_length));
        /** Exit if empty string or last element is NON_REF **/
        if(string_length == 0u) 
          break;
        else
        {
          if(insert_idx >= m_data.size())
            m_data.resize(2u*insert_idx+1u);
          char first_char = last_string[0];
          /** Exit if empty string or last element is NON_REF **/
          if(IS_NON_REF_ALLELE(first_char))
          {
            m_data[insert_idx++] = TILEDB_NON_REF_VARIANT_REPRESENTATION;
            break;
          }
          else
            m_data[insert_idx++] = std::move(last_string);
          /* DO NOT USE last_string beyond this point*/
          offset_idx += string_length + 1u; //the +1 is for going beyond the end '\0' 
        }
      }
      m_data.resize(insert_idx);
    }
    virtual std::vector<std::string>& get() { return m_data; }
    virtual const std::vector<std::string>& get() const { return m_data; }
    virtual void print(std::ostream& fptr) const
    {
      fptr << "[ ";
      for(auto& val : m_data)
        fptr << val << ",";
      fptr << "]";
    }
    virtual std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated)
    {
      size = m_data.size();
      char** data = new char*[size];
      for(auto i=0u;i<size;++i)
        data[i] = const_cast<char*>(m_data[i].c_str());
      *(reinterpret_cast<char***>(ptr)) = data;
      allocated = true;
      return std::type_index(typeid(char));
    }
    virtual VariantFieldBase* create_copy() const { return new VariantFieldALTData(*this); }
  private:
    std::vector<std::string> m_data;
};

/*
 * Base class of creator used by the VariantFieldFactory
 */
class VariantFieldCreatorBase {
  public:
    virtual std::unique_ptr<VariantFieldBase> Create() const = 0;
};

/*
 * Creator class that actually creates objects to hold field data. Invoked from within the factory
 */
template< class VariantFieldTy >
class VariantFieldCreator : public VariantFieldCreatorBase {
  public:
    VariantFieldCreator() {}
    virtual std::unique_ptr<VariantFieldBase> Create() const { return std::unique_ptr<VariantFieldBase>(new VariantFieldTy()); }
};

/**
 * Factory class to create VariantField objects given schema idx
 */
class VariantFieldFactory
{
  public:
    VariantFieldFactory() { clear(); }
    void clear() { m_schema_idx_to_creator.clear(); }
    void resize(unsigned num_fields_in_schema)
    {  m_schema_idx_to_creator.resize(num_fields_in_schema); }
    void Register(unsigned schema_idx, const std::shared_ptr<VariantFieldCreatorBase>& creator)
    {
      assert(schema_idx < m_schema_idx_to_creator.size());
      m_schema_idx_to_creator[schema_idx] = creator;
    }
    std::unique_ptr<VariantFieldBase> Create(unsigned schema_idx) const
    {
      assert(schema_idx < m_schema_idx_to_creator.size());
      assert(m_schema_idx_to_creator[schema_idx].get() != 0);
      return m_schema_idx_to_creator[schema_idx]->Create();
    }
  private:
    std::vector<std::shared_ptr<VariantFieldCreatorBase>> m_schema_idx_to_creator;
};
#endif
