#ifndef VARIANT_FIELD_DATA_H
#define VARIANT_FIELD_DATA_H

#include "query_field_data.h"
#include <memory>

/*
 * Base class for variant field data - not sure whether I will add any functionality here
 */
class BaseVariantFieldData : public QueryFieldData
{
  public:
    BaseVariantFieldData() 
      : QueryFieldData()
    { ; }
    virtual void copy_data_from_tile(const Tile&  base_tile, uint64_t element_idx, uint64_t num_elements) = 0;
    virtual void clear() { ; }
};
/*
 * Class that holds single element data 
 */
template<class DataType, class AttributeTileTy=const AttributeTile<DataType>>
class VariantFieldData : public BaseVariantFieldData
{
  public:
    VariantFieldData()
      : BaseVariantFieldData()
    { ; }
    virtual void copy_data_from_tile(const Tile& base_tile, uint64_t element_idx, uint64_t num_elements=0)
    {
      const AttributeTileTy& tile = static_cast<const AttributeTileTy&>(base_tile);
      tile.copy_payload(&m_data, element_idx, 1);
    }
    virtual const DataType& get() const { return m_data; }
  private:
    DataType m_data;
};
/*
 * Specialize class for string, AttributeTile type is <char> 
 */
template<>
class VariantFieldData<std::string, const AttributeTile<char>> : public BaseVariantFieldData
{
  public:
    VariantFieldData()
      : BaseVariantFieldData()
    { ; }
    virtual void clear() { m_data.clear(); }
    virtual void copy_data_from_tile(const Tile& base_tile, uint64_t element_idx, uint64_t num_elements=0)
    {
      const AttributeTile<char>& tile = static_cast<const AttributeTile<char>&>(base_tile);
      /** Assuming the string in the payload is null terminated, else dead.
       * This still copies the data into m_data, effectively strcpy, which hopefully is faster than any
       * custom code we write **/
      m_data = std::string(tile.get_payload_ptr(element_idx));
    }
    virtual const std::string& get() const { return m_data; }
  private:
    std::string m_data;
};
/*
 * Sub-class that holds vector data of basic types - int,float etc
 */
template<class DataType, class AttributeTileTy=const AttributeTile<DataType>>
class VariantFieldPrimitiveVectorData : public BaseVariantFieldData
{
  public:
    VariantFieldPrimitiveVectorData()
      : BaseVariantFieldData()
    { clear(); }
    virtual void clear() { m_data.clear(); }
    virtual void copy_data_from_tile(const Tile& base_tile, uint64_t element_idx, uint64_t num_elements)
    {
      const AttributeTileTy& tile = static_cast<const AttributeTileTy&>(base_tile);
      m_data.resize(num_elements);
      tile.copy_payload(&(m_data[0]), element_idx, num_elements);
    }
    virtual const std::vector<DataType>& get() const { return m_data; }
  private:
    std::vector<DataType> m_data;
};
/*
 * Special class for ALT field. ALT field is parsed in a weird way, hence treated in a special manner
 */
class VariantFieldALTData : public BaseVariantFieldData
{
  public:
    VariantFieldALTData()
      : BaseVariantFieldData()
    { ; }
    virtual void clear()
    {
      for(auto& s : m_data)
        s.clear();
      m_data.clear();
    }
    virtual void copy_data_from_tile(const Tile& base_tile, uint64_t element_idx, uint64_t num_elements=0)
    {
      m_data.clear();
      const AttributeTile<char>& tile = static_cast<const AttributeTile<char>&>(base_tile);
      uint64_t offset_idx = element_idx;
      /** Exit if empty string or last element is NON_REF **/
      while(1)
      {
        /* Assuming the string in the payload is null terminated, else dead.
         * This still copies the data into m_data, effectively calling strcpy which, hopefully, is faster than any
         * custom code we write **/
        std::string last_string = std::string(tile.get_payload_ptr(offset_idx));
        auto string_length = last_string.length();
        /** Exit if empty string or last element is NON_REF **/
        if(string_length == 0u) 
          break;
        else
        {
          char first_char = last_string[0];
          m_data.push_back(std::move(last_string));
          /* DO NOT USE last_string beyond this point*/
          /** Exit if empty string or last element is NON_REF **/
          if(IS_NON_REF_ALLELE(first_char))
            break;
          offset_idx += string_length + 1u; //the +1 is for going beyond the end '\0' 
        }
      }
    }
    virtual const std::vector<std::string>& get() const { return m_data; }
  private:
    std::vector<std::string> m_data;
};

/*
 * Base class of creator used by the VariantFieldFactory
 */
class VariantFieldCreatorBase {
  public:
    virtual std::unique_ptr<BaseVariantFieldData> Create() = 0;
};

/*
 * Creator class that actually creates objects to hold field data. Invoked from within the factory
 */
template< class VariantFieldDataTy >
class VariantFieldCreator : public VariantFieldCreatorBase {
  public:
    VariantFieldCreator() {}
    virtual std::unique_ptr<BaseVariantFieldData> Create() { return std::unique_ptr<BaseVariantFieldData>(new VariantFieldDataTy()); }
};

/**
 * Factory class to create VariantField objects given schema idx
 */
class VariantFieldFactory
{
  public:
    VariantFieldFactory() { m_schema_idx_to_creator.clear(); }
    void resize(unsigned num_fields_in_schema)
    {  m_schema_idx_to_creator.resize(num_fields_in_schema); }
    void Register(unsigned schema_idx, const std::shared_ptr<VariantFieldCreatorBase>& creator)
    {
      assert(schema_idx < m_schema_idx_to_creator.size());
      m_schema_idx_to_creator[schema_idx] = creator;
    }
    std::unique_ptr<BaseVariantFieldData> Create(unsigned schema_idx)
    {
      assert(schema_idx < m_schema_idx_to_creator.size());
      assert(m_schema_idx_to_creator[schema_idx].get() != 0);
      return m_schema_idx_to_creator[schema_idx]->Create();
    }
  private:
    std::vector<std::shared_ptr<VariantFieldCreatorBase>> m_schema_idx_to_creator;
};
#endif
