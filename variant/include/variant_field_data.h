#ifndef VARIANT_FIELD_DATA_H
#define VARIANT_FIELD_DATA_H

#include <memory>
#include <unordered_map>
#include "query_field_data.h"
#include "gt_common.h"
#include "cell_const_attr_iterator.h"
#include "special_values.h"

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
  VARIANT_FIELD_CHAR
};
extern std::unordered_map<std::type_index, VariantFieldTypeEnum> g_variant_field_type_index_to_enum;

template<class T>
bool is_tiledb_missing_value(const T value);

//For comparison
typedef union
{
  unsigned i;
  float f;
}fi_union;

typedef union
{
  uint64_t i;
  double d;
}di_union;

/*
 * Base class for variant field data - not sure whether I will add any functionality here
 */
class VariantFieldBase : public QueryFieldData
{
  public:
    VariantFieldBase() 
      : QueryFieldData()
    {
      m_subclass_type = VARIANT_FIELD_BASE;
      m_valid = false;
    }
    virtual ~VariantFieldBase() = default;
    virtual void copy_data_from_tile(const CellConstAttrIterator&  attr_iter, unsigned length_descriptor, unsigned num_elements) = 0;
    virtual void clear() { ; }
    virtual void print(std::ostream& fptr) const { ; }
    /* Get pointer(s) to data with number of elements */
    virtual std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated) = 0;
    /* Return type of data */
    virtual std::type_index get_element_type() const = 0;
    /* Create copy and return pointer - avoid using as much as possible*/
    virtual VariantFieldBase* create_copy() const = 0;
    //Resize this field - default do nothing
    virtual void resize(unsigned new_size) { ; }
    /* Return address of the offset-th element */
    virtual void* get_address(unsigned offset) { return 0; }
    /* Validity of field */
    void set_valid(bool value) { m_valid = value; }
    bool is_valid() const { return m_valid; }
  protected:
    enum VariantFieldTypesEnum
    {
      VARIANT_FIELD_BASE=0u,
      VARIANT_FIELD_DATA,
      VARIANT_FIELD_STRING,
      VARIANT_FIELD_PRIMITIVE_VECTOR,
      VARIANT_FIELD_ALT,
      NUM_VARIANT_FIELD_TYPES
    };
    unsigned m_subclass_type;   //enum from above
  private:
    bool m_valid;
};
/*
 * Class that holds single element data 
 */
template<class DataType>
class VariantFieldData : public VariantFieldBase
{
  public:
    VariantFieldData()
      : VariantFieldBase()
    { m_subclass_type = VARIANT_FIELD_DATA; }
    virtual ~VariantFieldData() = default;
    virtual void copy_data_from_tile(const CellConstAttrIterator&  attr_iter, uint64_t length_descriptor, unsigned num_elements)
    {
      auto base_ptr = (*attr_iter).operator const char*(); //const char*
      auto ptr = (*attr_iter).operator const DataType*(); //const DataType* ptr
      if(length_descriptor != BCF_VL_FIXED)     //variable length field, first 4 bytes are the length
        ptr = reinterpret_cast<const DataType*>(base_ptr + sizeof(int));
      m_data = *ptr;
    }
    virtual DataType& get() { return m_data; }
    virtual const DataType& get() const { return m_data; }
    virtual std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated)
    {
      size = 1u;
      *(reinterpret_cast<DataType**>(ptr)) = &m_data;
      allocated = false;
      return get_element_type();
    }
    virtual std::type_index get_element_type() const { return std::type_index(typeid(DataType)); }
    virtual VariantFieldBase* create_copy() const { return new VariantFieldData<DataType>(*this); }
    /* Return address of the offset-th element */
    virtual void* get_address(unsigned offset) { return reinterpret_cast<void*>(&m_data); }
  private:
    DataType m_data;
};
/*
 * Specialize class for string, AttributeTile type is <char> 
 */
template<>
class VariantFieldData<std::string> : public VariantFieldBase
{
  public:
    VariantFieldData()
      : VariantFieldBase()
    { m_subclass_type = VARIANT_FIELD_STRING; }
    virtual ~VariantFieldData() = default;
    virtual void clear() { m_data.clear(); }
    virtual void copy_data_from_tile(const CellConstAttrIterator&  attr_iter, unsigned length_descriptor, unsigned num_elements)
    {
      auto base_ptr = (*attr_iter).operator const char*(); //const char*
      auto ptr = (*attr_iter).operator const char*(); //const char* pointer
      if(length_descriptor != BCF_VL_FIXED)     //variable length field, first 4 bytes are the length
      {
        auto num_elements_ptr = reinterpret_cast<const int*>(base_ptr);
        num_elements = *num_elements_ptr;
        ptr = static_cast<const char*>(base_ptr + sizeof(int));
      }
      m_data = std::move(std::string(ptr, num_elements));
      bool is_missing_flag = true;
      for(auto val : m_data)
        if(val != NULL_CHAR)
        {
          is_missing_flag = false;
          break;
        }
      //Whole field is missing, clear and invalidate
      //Note that 0-length fields are invalidated based on the logic in this function
      if(is_missing_flag)
      {
        set_valid(false);
        m_data.clear();
      }
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
      return get_element_type();
    }
    virtual std::type_index get_element_type() const { return std::type_index(typeid(char)); }
    virtual VariantFieldBase* create_copy() const { return new VariantFieldData<std::string>(*this); }
    /* Return address of the offset-th element */
    virtual void* get_address(unsigned offset) { return reinterpret_cast<void*>(&m_data); }
  private:
    std::string m_data;
};
//Assigned name for string type
typedef VariantFieldData<std::string> VariantFieldString;
/*
 * Sub-class that holds vector data of basic types - int,float etc
 */
template<class DataType>
class VariantFieldPrimitiveVectorData : public VariantFieldBase
{
  public:
    VariantFieldPrimitiveVectorData()
      : VariantFieldBase()
    { 
      m_subclass_type = VARIANT_FIELD_PRIMITIVE_VECTOR;
      clear(); 
    }
    virtual ~VariantFieldPrimitiveVectorData() = default;
    virtual void clear() { m_data.clear(); }
    virtual void copy_data_from_tile(const CellConstAttrIterator&  attr_iter, unsigned length_descriptor, unsigned num_elements)
    {
      auto base_ptr = (*attr_iter).operator const char*(); //const char*
      auto ptr = (*attr_iter).operator const DataType*(); //const DataType* ptr
      if(length_descriptor != BCF_VL_FIXED)     //variable length field, first 4 bytes are the length
      {
        auto num_elements_ptr = reinterpret_cast<const int*>(base_ptr);
        num_elements = *num_elements_ptr;
        ptr = reinterpret_cast<const DataType*>(base_ptr + sizeof(int));
      }
      m_data.resize(num_elements);
      memcpy(&(m_data[0]), ptr, num_elements*sizeof(DataType));
      bool is_missing_flag = true;
      for(auto val : m_data)
        if(!is_tiledb_missing_value<DataType>(val))
        {
          is_missing_flag = false;
          break;
        }
      //Whole field is missing, clear and invalidate
      //Note that 0-length fields are invalidated based on the logic in this function
      if(is_missing_flag)
      {
        set_valid(false);
        m_data.clear();
      }
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
      return get_element_type();
    }
    virtual std::type_index get_element_type() const { return std::type_index(typeid(DataType)); }
    virtual VariantFieldBase* create_copy() const { return new VariantFieldPrimitiveVectorData<DataType>(*this); }
    virtual void resize(unsigned new_size) { m_data.resize(new_size); }
    /* Return address of the offset-th element */
    virtual void* get_address(unsigned offset)
    {
      assert(offset < m_data.size());
      return reinterpret_cast<void*>(&(m_data[offset]));
    }
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
    { 
      m_subclass_type = VARIANT_FIELD_ALT;
      clear();
    }
    virtual ~VariantFieldALTData() = default;
    virtual void clear()
    {
      for(auto& s : m_data)
        s.clear();
      m_data.clear();
    }
    virtual void copy_data_from_tile(const CellConstAttrIterator&  attr_iter, unsigned length_descriptor, unsigned num_elements)
    {
      auto base_ptr = (*attr_iter).operator const char*(); //const char*
      auto ptr = (*attr_iter).operator const char*(); //const char* pointer
      if(length_descriptor != BCF_VL_FIXED)     //variable length field, first 4 bytes are the length
      {
        auto num_elements_ptr = reinterpret_cast<const int*>(base_ptr);
        num_elements = *num_elements_ptr;
        ptr = static_cast<const char*>(base_ptr + sizeof(int));
      }
      //Create copy for use in strtok
      char* tmp = strndup(ptr, num_elements);
      assert(strnlen(tmp, num_elements+1) == num_elements); //last element should be 0
      //Tokenize
      char* saveptr = 0;
      char* argptr = tmp;
      m_data.clear();
      while(auto curr_token = strtok_r(argptr, TILEDB_ALT_ALLELE_SEPARATOR, &saveptr))
      {
        m_data.push_back(curr_token);
        argptr = 0;
      }
      free(tmp);
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
    virtual std::type_index get_element_type() const { return std::type_index(typeid(std::string)); }   //each element is a string
    virtual VariantFieldBase* create_copy() const { return new VariantFieldALTData(*this); }
    virtual void resize(unsigned new_size) { m_data.resize(new_size); }
    /* Return address of the offset-th element */
    virtual void* get_address(unsigned offset)
    {
      assert(offset < m_data.size());
      return reinterpret_cast<void*>(&(m_data[offset]));
    }
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
