#ifndef VARIANT_FIELD_DATA_H
#define VARIANT_FIELD_DATA_H

#include <memory>
#include "headers.h"
#include "gt_common.h"
#include "variant_cell.h"

class UnknownAttributeTypeException : public std::exception {
  public:
    UnknownAttributeTypeException(const std::string m="Unknown type of queried attribute") : msg_(m) { ; }
    ~UnknownAttributeTypeException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class VariantFieldTypeUtil
{
  public:
    static size_t size(VariantFieldTypeEnum enum_val);
    static size_t size(const std::type_index& type_index)
    {
      auto iter = g_variant_field_type_index_to_enum.find(type_index);
      if(iter == g_variant_field_type_index_to_enum.end())
        throw UnknownAttributeTypeException(std::string("Unhandled attribute type ")+type_index.name());
      return size((*iter).second);
    }
    static int get_tiledb_type_for_variant_field_type(const std::type_index& type_index)
    {
      auto iter = g_variant_field_type_index_to_tiledb_type.find(type_index);
      if(iter == g_variant_field_type_index_to_tiledb_type.end())
       throw UnknownAttributeTypeException(std::string("No TileDB type found for attribute ")+type_index.name());
      return (*iter).second;
    }
    static std::type_index get_variant_field_type_for_tiledb_type(const int tiledb_type)
    {
      assert(static_cast<size_t>(tiledb_type) < g_tiledb_type_to_variant_field_type_index.size());
      auto type_idx = g_tiledb_type_to_variant_field_type_index[tiledb_type];
      assert(g_variant_field_type_index_to_tiledb_type.find(type_idx) != g_variant_field_type_index_to_tiledb_type.end()
          && g_variant_field_type_index_to_tiledb_type[type_idx] == tiledb_type);
      return type_idx;
    }
};

template<class T>
bool is_tiledb_missing_value(const T value);

#define RESIZE_BINARY_SERIALIZATION_BUFFER_IF_NEEDED(buffer, offset, add_size) \
      if(offset + add_size > buffer.size()) \
        buffer.resize(offset + add_size + 1024u);
/*
 * Base class for variant field data - not sure whether I will add any functionality here
 */
class VariantFieldBase
{
  public:
    VariantFieldBase() 
    {
      m_subclass_type = VARIANT_FIELD_BASE;
      m_valid = false;
    }
    virtual ~VariantFieldBase() = default;
    virtual void copy_data_from_tile(const BufferVariantCell::FieldsIter&  attr_iter) = 0;
    virtual void clear() { ; }
    virtual void print(std::ostream& fptr) const { ; }
    virtual void print_Cotton_JSON(std::ostream& fptr) const { ; }
    virtual void binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const = 0;
    virtual void binary_deserialize(const char* buffer, uint64_t& offset, unsigned length_descriptor, unsigned num_elements) = 0;
    /* Get pointer(s) to data with number of elements */
    virtual std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated) = 0;
    /* Get raw pointer(s) to data */
    virtual const void* get_raw_pointer() const = 0;
    /* Return type of data */
    virtual std::type_index get_element_type() const = 0;
    /* Return #elements */
    virtual size_t length() const = 0;
    /* Create copy and return pointer - avoid using as much as possible*/
    virtual VariantFieldBase* create_copy() const = 0;
    /* Copy from src ptr*/
    virtual void copy_from(const VariantFieldBase* base_src)
    {
      m_valid = base_src->m_valid;
      m_subclass_type = base_src->m_subclass_type;
    }
    //Resize this field - default do nothing
    virtual void resize(unsigned new_size) { ; }
    /* Return address of the offset-th element */
    virtual void* get_address(unsigned offset) { return 0; }
    /* Validity of field */
    void set_valid(bool value) { m_valid = value; }
    bool is_valid() const { return m_valid; }
  protected:
    enum VariantFieldClassTypeEnum
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
    virtual void copy_data_from_tile(const BufferVariantCell::FieldsIter&  attr_iter)
    {
      auto base_ptr = attr_iter.operator*<char>(); //const char*
      uint64_t offset = 0ull;
      //Set length descriptor to BCF_VL_FIXED as attr_iter will return pointer to data directly
      //attr_iter will have consumed the length field internally
      binary_deserialize(base_ptr, offset, BCF_VL_FIXED, attr_iter.get_field_length());
    }
    virtual void binary_deserialize(const char* buffer, uint64_t& offset, unsigned length_descriptor, unsigned num_elements)
    {
      auto base_ptr = buffer + offset; //const char*
      auto ptr = reinterpret_cast<const DataType*>(base_ptr); //const DataType* ptr
      if(length_descriptor != BCF_VL_FIXED)     //variable length field, first 4 bytes are the length
      {
        ptr = reinterpret_cast<const DataType*>(base_ptr + sizeof(int));
        offset += sizeof(int);
      }
      m_data = *ptr;
      offset += sizeof(DataType);
    }
    virtual void binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const
    {
      uint64_t add_size = sizeof(DataType);
      RESIZE_BINARY_SERIALIZATION_BUFFER_IF_NEEDED(buffer, offset, add_size);
      *(reinterpret_cast<DataType*>(&(buffer[offset]))) = m_data;
      offset += sizeof(DataType);
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
    virtual const void* get_raw_pointer() const  { return reinterpret_cast<void*>(&m_data); }
    virtual std::type_index get_element_type() const { return std::type_index(typeid(DataType)); }
    virtual size_t length() const { return 1u; };
    virtual VariantFieldBase* create_copy() const { return new VariantFieldData<DataType>(*this); }
    virtual void copy_from(const VariantFieldBase* base_src)
    {
      VariantFieldBase::copy_from(base_src);
      auto src = dynamic_cast<const VariantFieldData<DataType>*>(base_src);
      assert(src);
      m_data = src->m_data;
    }
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
    virtual void copy_data_from_tile(const BufferVariantCell::FieldsIter&  attr_iter)
    {
      auto base_ptr = attr_iter.operator*<char>(); //const char*
      uint64_t offset = 0ull;
      //Set length descriptor to BCF_VL_FIXED as attr_iter will return pointer to data directly
      //attr_iter will have consumed the length field internally
      binary_deserialize(base_ptr, offset, BCF_VL_FIXED, attr_iter.get_field_length());
    }
    virtual void binary_deserialize(const char* buffer, uint64_t& offset, unsigned length_descriptor, unsigned num_elements)
    {
      auto base_ptr = buffer + offset; //const char* pointer
      auto ptr = base_ptr;      //const char* pointer
      if(length_descriptor != BCF_VL_FIXED)     //variable length field, first 4 bytes are the length
      {
        auto num_elements_ptr = reinterpret_cast<const int*>(ptr);
        num_elements = *num_elements_ptr;
        ptr = static_cast<const char*>(base_ptr + sizeof(int));
        offset += sizeof(int);
      }
      m_data.resize(num_elements);
      memcpy(&(m_data[0]), ptr, num_elements*sizeof(char));
      bool is_missing_flag = true;
      for(auto val : m_data)
        if(!is_tiledb_missing_value<char>(val))
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
      offset += num_elements*sizeof(char);
    }
    virtual std::string& get()  { return m_data; }
    virtual const std::string& get() const { return m_data; }
    virtual void print(std::ostream& fptr) const { fptr << m_data; }
    virtual void print_Cotton_JSON(std::ostream& fptr) const { fptr << "\"" << m_data << "\"" ; }
    virtual void binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const
    {
      //string length + contents
      unsigned str_length = m_data.length();
      uint64_t add_size = sizeof(int) + str_length;
      RESIZE_BINARY_SERIALIZATION_BUFFER_IF_NEEDED(buffer, offset, add_size);
      //string length
      *(reinterpret_cast<int*>(&(buffer[offset]))) = str_length;
      offset += sizeof(int);
      //string contents
      memcpy(&(buffer[offset]), &(m_data[0]), str_length);
      offset += str_length;
    }
    virtual std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated)
    {
      size = 1u;
      char** data = new char*;
      data[0] = const_cast<char*>(m_data.c_str());
      *(reinterpret_cast<char***>(ptr)) = data;
      allocated = true;
      return get_element_type();
    }
    virtual const void* get_raw_pointer() const  { return reinterpret_cast<const void*>(m_data.c_str()); }
    virtual std::type_index get_element_type() const { return std::type_index(typeid(char)); }
    virtual size_t length() const { return m_data.length(); };
    virtual VariantFieldBase* create_copy() const { return new VariantFieldData<std::string>(*this); }
    virtual void copy_from(const VariantFieldBase* base_src)
    {
      VariantFieldBase::copy_from(base_src);
      auto src = dynamic_cast<const VariantFieldData<std::string>*>(base_src);
      assert(src);
      m_data.resize(src->m_data.size());
      if(m_data.size())
        memcpy(&(m_data[0]), &(src->m_data[0]), m_data.size()*sizeof(char));
    }
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
    virtual void copy_data_from_tile(const BufferVariantCell::FieldsIter&  attr_iter)
    {
      auto base_ptr = attr_iter.operator*<char>(); //const char*
      uint64_t offset = 0ull;
      //Set length descriptor to BCF_VL_FIXED as attr_iter will return pointer to data directly
      //attr_iter will have consumed the length field internally
      binary_deserialize(base_ptr, offset, BCF_VL_FIXED, attr_iter.get_field_length());
      m_length_descriptor = attr_iter.is_variable_length_field() ? BCF_VL_VAR : BCF_VL_FIXED;
    }
    virtual void binary_deserialize(const char* buffer, uint64_t& offset, unsigned length_descriptor, unsigned num_elements)
    {
      auto base_ptr = buffer + offset; //const char*
      auto ptr = reinterpret_cast<const DataType*>(base_ptr); //const DataType* ptr
      m_length_descriptor = length_descriptor;
      if(length_descriptor != BCF_VL_FIXED)     //variable length field, first 4 bytes are the length
      {
        auto num_elements_ptr = reinterpret_cast<const int*>(base_ptr);
        num_elements = *num_elements_ptr;
        ptr = reinterpret_cast<const DataType*>(base_ptr + sizeof(int));
        offset += sizeof(int);
      }
      m_data.resize(num_elements);
      unsigned data_size = num_elements*sizeof(DataType);
      memcpy(&(m_data[0]), ptr, data_size);
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
      offset += data_size; 
    }
    virtual std::vector<DataType>& get()  { return m_data; }
    virtual const std::vector<DataType>& get() const { return m_data; }
    virtual void print(std::ostream& fptr) const
    {
      fptr << "[ ";
      auto first_elem = true;
      for(auto val : m_data)
      {
        if(first_elem)
        {
          fptr << val;
          first_elem = false;
        }
        else
          fptr << "," << val;
      }
      fptr << " ]";
    }
    virtual void print_Cotton_JSON(std::ostream& fptr) const
    {
      //Variable length field or #elements > 1, print JSON list
      if(m_length_descriptor != BCF_VL_FIXED || m_data.size() > 1u)
        print(fptr);
      else      //single element field
        if(m_data.size() > 0u)
          fptr << m_data[0];
        else
          fptr << "null";
    }
    virtual void binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const
    {
      //Data contents
      unsigned data_length = m_data.size()*sizeof(DataType);
      //Add length field, if var sized field
      uint64_t add_size = ((m_length_descriptor == BCF_VL_FIXED) ? 0u : sizeof(int)) + data_length;
      RESIZE_BINARY_SERIALIZATION_BUFFER_IF_NEEDED(buffer, offset, add_size);
      //Num elements
      if(m_length_descriptor != BCF_VL_FIXED)
      {
        *(reinterpret_cast<int*>(&(buffer[offset]))) = m_data.size();
        offset += sizeof(int);
      }
      //data contents
      memcpy(&(buffer[offset]), &(m_data[0]), data_length);
      offset += data_length;
    }
    virtual std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated)
    {
      size = m_data.size();
      *(reinterpret_cast<DataType**>(ptr)) = (size > 0) ? &(m_data[0]) : nullptr;
      allocated = false;
      return get_element_type();
    }
    virtual const void* get_raw_pointer() const  { return reinterpret_cast<const void*>(m_data.size() ? &(m_data[0]) : 0); }
    virtual std::type_index get_element_type() const { return std::type_index(typeid(DataType)); }
    virtual size_t length() const { return m_data.size(); }
    virtual VariantFieldBase* create_copy() const { return new VariantFieldPrimitiveVectorData<DataType>(*this); }
    virtual void copy_from(const VariantFieldBase* base_src)
    {
      VariantFieldBase::copy_from(base_src);
      auto src = dynamic_cast<const VariantFieldPrimitiveVectorData<DataType>*>(base_src);
      assert(src);
      m_length_descriptor = src->m_length_descriptor;
      m_data.resize(src->m_data.size());
      if(m_data.size())
        memcpy(&(m_data[0]), &(src->m_data[0]), m_data.size()*sizeof(DataType));
    }
    virtual void resize(unsigned new_size) { m_data.resize(new_size); }
    /* Return address of the offset-th element */
    virtual void* get_address(unsigned offset)
    {
      assert(offset < m_data.size());
      return reinterpret_cast<void*>(&(m_data[offset]));
    }
  private:
    std::vector<DataType> m_data;
    unsigned m_length_descriptor;
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
    virtual void copy_data_from_tile(const BufferVariantCell::FieldsIter&  attr_iter)
    {
      auto base_ptr = attr_iter.operator*<char>(); //const char*
      uint64_t offset = 0ull;
      //Set length descriptor to BCF_VL_FIXED as attr_iter will return pointer to data directly
      //attr_iter will have consumed the length field internally
      binary_deserialize(base_ptr, offset, BCF_VL_FIXED, attr_iter.get_field_length());
    }
    virtual void binary_deserialize(const char* buffer, uint64_t& offset, unsigned length_descriptor, unsigned num_elements)
    {
      auto base_ptr = buffer + offset; //const char*
      auto ptr = base_ptr; //const char* pointer
      if(length_descriptor != BCF_VL_FIXED)     //variable length field, first 4 bytes are the length
      {
        auto num_elements_ptr = reinterpret_cast<const int*>(base_ptr);
        num_elements = *num_elements_ptr;
        ptr = static_cast<const char*>(base_ptr + sizeof(int));
        offset += sizeof(int);
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
      offset += num_elements*sizeof(char);
    }
    virtual std::vector<std::string>& get() { return m_data; }
    virtual const std::vector<std::string>& get() const { return m_data; }
    virtual void print(std::ostream& fptr) const
    {
      fptr << "[ ";
      auto first_elem = true;
      for(auto& val : m_data)
      {
        auto& ptr = IS_NON_REF_ALLELE(val) ? g_vcf_NON_REF : val;
        if(first_elem)
        {
          fptr << "\"" << ptr << "\"";
          first_elem = false;
        }
        else
          fptr << ",\"" << ptr << "\"";
      }
      fptr << " ]";
    }
    virtual void print_Cotton_JSON(std::ostream& fptr) const { print(fptr); }
    virtual void binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const
    {
      //string length
      uint64_t add_size = sizeof(int);
      RESIZE_BINARY_SERIALIZATION_BUFFER_IF_NEEDED(buffer, offset, add_size);
      //store location at which string length must be stored
      auto str_length_offset = offset;
      offset += sizeof(int);
      //location at which ALT string begins
      auto str_begin_offset = offset;
      bool first_elem = true;
      for(auto& val : m_data)
      {
        auto str_size = val.length();
        if(first_elem)
        {
          add_size = str_size;
          RESIZE_BINARY_SERIALIZATION_BUFFER_IF_NEEDED(buffer, offset, add_size);
          first_elem = false;
        }
        else    //add '|' separator
        {
          add_size = str_size + sizeof(char); //+1 for ALT separator
          RESIZE_BINARY_SERIALIZATION_BUFFER_IF_NEEDED(buffer, offset, add_size);
          *(reinterpret_cast<char*>(&(buffer[offset]))) = *TILEDB_ALT_ALLELE_SEPARATOR;
          offset += sizeof(char);
        }
        memcpy(&(buffer[offset]), val.c_str(), str_size);
        offset += str_size;
      }
      //string length
      *(reinterpret_cast<int*>(&(buffer[str_length_offset]))) = offset - str_begin_offset;
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
    virtual const void* get_raw_pointer() const  { return reinterpret_cast<const void*>(m_data.size() ? &(m_data[0]) : 0); }
    virtual std::type_index get_element_type() const { return std::type_index(typeid(std::string)); }   //each element is a string
    virtual size_t length() const { return m_data.size(); }
    virtual VariantFieldBase* create_copy() const { return new VariantFieldALTData(*this); }
    virtual void copy_from(const VariantFieldBase* base_src)
    {
      VariantFieldBase::copy_from(base_src);
      auto src = dynamic_cast<const VariantFieldALTData*>(base_src);
      assert(src);
      m_data.resize(src->m_data.size());
      for(auto i=0u;i<m_data.size();++i)
      {
        auto& curr_dst = m_data[i];
        auto& curr_src = src->m_data[i];
        curr_dst.resize(curr_src.size());
        memcpy(&(curr_dst[0]), &(curr_src[0]), curr_src.size()*sizeof(char));
      }
    }
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
