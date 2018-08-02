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
    static size_t size(const int bcf_ht_type);
    static size_t size(const std::type_index& type_index)
    {
      auto iter = g_variant_field_type_index_to_bcf_ht_type.find(type_index);
      if(iter == g_variant_field_type_index_to_bcf_ht_type.end())
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
   
    static int get_bcf_ht_type_for_variant_field_type(const std::type_index& type_index)
    {
      auto iter = g_variant_field_type_index_to_bcf_ht_type.find(type_index);
      if(iter == g_variant_field_type_index_to_bcf_ht_type.end())
        throw UnknownAttributeTypeException(std::string("Unhandled attribute type ")+type_index.name());
      return (*iter).second;
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
    VariantFieldBase(const bool is_variable_length_field) 
    {
      m_subclass_type = VARIANT_FIELD_BASE;
      m_valid = false;
      m_is_variable_length_field = is_variable_length_field;
      m_cell_idx = 0;
    }
    virtual ~VariantFieldBase() = default;
    void copy_data_from_tile(const BufferVariantCell::FieldsIter& attr_iter)
    {
      auto base_ptr = attr_iter.operator*<char>(); //const char*
      uint64_t offset = 0ull;
      //buffer does not have length serialized
      binary_deserialize(base_ptr, offset, false, attr_iter.get_field_length());
    }
    virtual void print(std::ostream& fptr) const  = 0;
    virtual void print_csv(std::ostream& fptr) const  = 0;
    virtual void print_Cotton_JSON(std::ostream& fptr) const { ; }
    virtual void binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const = 0; 
    virtual void binary_deserialize(const char* buffer, uint64_t& offset,
        const bool is_length_in_buffer, unsigned num_elements) = 0;
    /* Get pointer(s) to data with number of elements */
    virtual std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated) = 0;
    /* Get raw pointer(s) to data */
    virtual const void* get_raw_pointer() const = 0;
    /* Return #elements */
    virtual size_t length() const = 0;
    /* Create copy and return pointer - avoid using as much as possible*/
    virtual VariantFieldBase* create_copy() const = 0;
    /* Copy from src ptr*/
    virtual void copy_from(const VariantFieldBase* base_src)
    {
      m_valid = base_src->m_valid;
      m_is_variable_length_field = base_src->m_is_variable_length_field;
      m_subclass_type = base_src->m_subclass_type;
      m_cell_idx = base_src->m_cell_idx;
    }
    //Resize this field - default do nothing
    virtual void resize(unsigned new_size) { ; }
    /* Return address of the offset-th element */
    virtual void* get_address(unsigned offset) { return 0; }
    /* Validity of field */
    void set_valid(bool value) { m_valid = value; }
    bool is_valid() const { return m_valid; }
    /*
     * Columnar field info
     */
    void set_columnar_field_object(GenomicsDBColumnarField* columnar_field, const uint64_t cell_idx)
    {
      m_cell_idx = cell_idx;
    }
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
    bool m_is_variable_length_field;
    unsigned m_subclass_type;   //enum from above
    uint64_t m_cell_idx;
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
    VariantFieldData(const bool is_variable_length_field=false)
      : VariantFieldBase(false)
    { m_subclass_type = VARIANT_FIELD_DATA; }
    virtual ~VariantFieldData() = default;
    virtual void binary_deserialize(const char* buffer, uint64_t& offset,
        const bool is_length_in_buffer, unsigned num_elements)
    {
      auto base_ptr = buffer + offset; //const char*
      auto ptr = reinterpret_cast<const DataType*>(base_ptr); //const DataType* ptr
      if(is_length_in_buffer)     //length specified in buffer
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
    DataType& get() { return m_data; }
    const DataType& get() const { return m_data; }
    virtual std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated)
    {
      size = 1u;
      *(reinterpret_cast<DataType**>(ptr)) = &m_data;
      allocated = false;
      return get_element_type();
    }
    virtual const void* get_raw_pointer() const  { return reinterpret_cast<void*>(&m_data); }
    std::type_index get_element_type() const { return std::type_index(typeid(DataType)); }
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
    virtual void print(std::ostream& fptr) const { fptr << m_data; }
    virtual void print_csv(std::ostream& fptr) const { fptr << m_data; }
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
    VariantFieldData(const bool is_variable_length_field=true)
      : VariantFieldBase(true)
    { m_subclass_type = VARIANT_FIELD_STRING; }
    virtual ~VariantFieldData() = default;
    virtual void binary_deserialize(const char* buffer, uint64_t& offset,
        const bool is_length_in_buffer, unsigned num_elements)
    {
      auto base_ptr = buffer + offset; //const char* pointer
      auto ptr = base_ptr;      //const char* pointer
      if(is_length_in_buffer)     //length specified in buffer
      {
        auto num_elements_ptr = reinterpret_cast<const int*>(ptr);
        num_elements = *num_elements_ptr;
        ptr = static_cast<const char*>(base_ptr + sizeof(int));
        offset += sizeof(int);
      }
      m_data.resize(num_elements);
      memcpy_s(&(m_data[0]), num_elements*sizeof(char), ptr, num_elements*sizeof(char));
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
    virtual void print(std::ostream& fptr) const { fptr << "\"" << m_data << "\""; }
    virtual void print_csv(std::ostream& fptr) const { fptr << m_data; }
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
      memcpy_s(&(buffer[offset]), str_length, &(m_data[0]), str_length);
      offset += str_length;
    }
    std::string& get() { return m_data; }
    const std::string& get() const { return m_data; }
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
    std::type_index get_element_type() const { return std::type_index(typeid(char)); }
    virtual size_t length() const { return m_data.length(); };
    virtual VariantFieldBase* create_copy() const { return new VariantFieldData<std::string>(*this); }
    virtual void copy_from(const VariantFieldBase* base_src)
    {
      VariantFieldBase::copy_from(base_src);
      auto src = dynamic_cast<const VariantFieldData<std::string>*>(base_src);
      assert(src);
      m_data.resize(src->m_data.size());
      if(m_data.size())
        memcpy_s(&(m_data[0]), m_data.size()*sizeof(char), &(src->m_data[0]), m_data.size()*sizeof(char));
    }
    /* Return address of the offset-th element */
    virtual void* get_address(unsigned offset) { return reinterpret_cast<void*>(&m_data); }
  private:
    std::string m_data;
};
//Assigned name for string type
typedef VariantFieldData<std::string> VariantFieldString;

//For some common functions
class VariantFieldPrimitiveVectorDataBase : public VariantFieldBase
{
  public:
    VariantFieldPrimitiveVectorDataBase(const bool is_variable_length_field)
      : VariantFieldBase(is_variable_length_field)
    { 
      m_subclass_type = VARIANT_FIELD_PRIMITIVE_VECTOR;
    }
    virtual void binary_deserialize(const char* buffer, uint64_t& offset,
        const bool is_length_in_buffer, unsigned num_elements)
    {
      auto base_ptr = buffer + offset; //const char*
      auto ptr = base_ptr;
      if(is_length_in_buffer)     //length specified in buffer
      {
        auto num_elements_ptr = reinterpret_cast<const int*>(base_ptr);
        num_elements = *num_elements_ptr;
        ptr = base_ptr + sizeof(int);
        offset += sizeof(int);
      }
      copy_data_into_vector(ptr, num_elements);
      offset += num_elements*element_size();
    }
    virtual void binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const
    {
      //Data contents
      unsigned data_length = length()*element_size();
      //Add length field, if var sized field
      uint64_t add_size = (m_is_variable_length_field ? sizeof(int) : 0u) + data_length;
      RESIZE_BINARY_SERIALIZATION_BUFFER_IF_NEEDED(buffer, offset, add_size);
      //Num elements
      if(m_is_variable_length_field)
      {
        *(reinterpret_cast<int*>(&(buffer[offset]))) = length();
        offset += sizeof(int);
      }
      //data contents
      memcpy_s(&(buffer[offset]), data_length, get_raw_pointer(), data_length);
      offset += data_length;
    }
    virtual void copy_data_into_vector(const char* ptr, const size_t num_bytes) = 0;
    virtual size_t element_size() const = 0;
};

/*
 * Sub-class that holds vector data of basic types - int,float etc
 */
template<class DataType, class PrintType=DataType>
class VariantFieldPrimitiveVectorData : public VariantFieldPrimitiveVectorDataBase
{
  public:
    VariantFieldPrimitiveVectorData(const bool is_variable_length_field)
      : VariantFieldPrimitiveVectorDataBase(is_variable_length_field)
    { 
    }
    virtual ~VariantFieldPrimitiveVectorData() = default;
    std::type_index get_element_type() const { return std::type_index(typeid(DataType)); }
    virtual VariantFieldBase* create_copy() const { return new VariantFieldPrimitiveVectorData<DataType, PrintType>(*this); }
    size_t element_size() const { return sizeof(DataType); }
    void copy_data_into_vector(const char* buffer, const size_t num_elements)
    {
      m_data.resize(num_elements);
      unsigned data_size = num_elements*sizeof(DataType);
      memcpy_s(&(m_data[0]), data_size, buffer, data_size);
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
    std::vector<DataType>& get()  { return m_data; }
    const std::vector<DataType>& get() const { return m_data; }
    void print(std::ostream& fptr) const
    {
      fptr << "[ ";
      auto first_elem = true;
      for(auto val : m_data)
      {
        if(first_elem)
        {
          fptr << static_cast<PrintType>(val);
          first_elem = false;
        }
        else
          fptr << "," << static_cast<PrintType>(val);
      }
      fptr << " ]";
    }
    void print_phased_GT_field(std::ostream& fptr)
    {
      fptr << "[ ";
      auto first_elem = true;
      auto is_phase_element = false;
      for(auto val : m_data)
      {
        if(first_elem)
        {
          fptr << static_cast<PrintType>(val);
          first_elem = false;
        }
        else
        {
          fptr << ",";
          if(is_phase_element)
          {
            if(val)
              fptr << "|";
            else
              fptr << "\\";
          }
          else
            fptr << static_cast<PrintType>(val);
        }
        is_phase_element = !is_phase_element; //flip
      }
      fptr << " ]";
    }

    void print_csv(std::ostream& fptr) const
    {
      if(m_is_variable_length_field)
        fptr << m_data.size() << ",";
      //Print blanks for invalid cells
      auto first_elem = true;
      for(auto val : m_data)
      {
        if(first_elem)
        {
          fptr << static_cast<PrintType>(val);
          first_elem = false;
        }
        else
          fptr << "," << static_cast<PrintType>(val);
      }
    }
    void print_Cotton_JSON(std::ostream& fptr) const
    {
      //Variable length field or #elements > 1, print JSON list
      if(m_is_variable_length_field || m_data.size() > 1u)
        print(fptr);
      else      //single element field
        if(m_data.size() > 0u)
          fptr << m_data[0];
        else
          fptr << "null";
    } 
    std::type_index get_C_pointers(unsigned& size, void** ptr, bool& allocated)
    {
      size = m_data.size();
      *(reinterpret_cast<DataType**>(ptr)) = (size > 0) ? &(m_data[0]) : nullptr;
      allocated = false;
      return get_element_type();
    }
    const void* get_raw_pointer() const
    {
      return reinterpret_cast<const void*>(m_data.size() ? &(m_data[0]) : 0);
    }
    size_t length() const
    {
      return m_data.size();
    }
    void copy_from(const VariantFieldBase* base_src)
    {
      VariantFieldBase::copy_from(base_src);
      auto src = dynamic_cast<const VariantFieldPrimitiveVectorData<DataType, PrintType>*>(base_src);
      assert(src);
      m_data.resize(src->m_data.size());
      if(m_data.size())
        memcpy_s(&(m_data[0]), m_data.size()*sizeof(DataType), &(src->m_data[0]), m_data.size()*sizeof(DataType));
    }
    void resize(unsigned new_size)
    {
      m_data.resize(new_size);
    }
    /*Return address of the offset-th element*/
    void* get_address(unsigned offset)
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
    VariantFieldALTData(const bool is_variable_length_field=true)
      : VariantFieldBase(true)
    { 
      m_subclass_type = VARIANT_FIELD_ALT;
    }
    virtual ~VariantFieldALTData() = default;
    virtual void binary_deserialize(const char* buffer, uint64_t& offset,
        const bool is_length_in_buffer, unsigned num_elements)
    {
      auto base_ptr = buffer + offset; //const char*
      auto ptr = base_ptr; //const char* pointer
      if(is_length_in_buffer)     //length specified in buffer
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
    virtual void print_csv(std::ostream& fptr) const
    {
      auto first_elem = true;
      for(auto& val : m_data)
      {
        if(first_elem)
        {
          fptr << val;
          first_elem = false;
        }
        else
          fptr << "," << val;
      }
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
        memcpy_s(&(buffer[offset]), str_size, val.c_str(), str_size);
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
    std::type_index get_element_type() const { return std::type_index(typeid(std::string)); }   //each element is a string
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
        memcpy_s(&(curr_dst[0]), curr_src.size()*sizeof(char), &(curr_src[0]), curr_src.size()*sizeof(char));
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
    virtual std::unique_ptr<VariantFieldBase> Create(const bool is_variable_length_field) const = 0;
};

/*
 * Creator class that actually creates objects to hold field data. Invoked from within the factory
 */
template< class VariantFieldTy >
class VariantFieldCreator : public VariantFieldCreatorBase {
  public:
    VariantFieldCreator() {}
    virtual std::unique_ptr<VariantFieldBase> Create(const bool is_variable_length_field) const
    {
      return std::unique_ptr<VariantFieldBase>(new VariantFieldTy(is_variable_length_field));
    }
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
    void Register(const unsigned schema_idx, const std::shared_ptr<VariantFieldCreatorBase>& creator)
    {
      assert(schema_idx < m_schema_idx_to_creator.size());
      m_schema_idx_to_creator[schema_idx] = creator;
    }
    std::unique_ptr<VariantFieldBase> Create(unsigned schema_idx, const bool is_variable_length_field) const
    {
      assert(schema_idx < m_schema_idx_to_creator.size());
      assert(m_schema_idx_to_creator[schema_idx].get() != 0);
      return m_schema_idx_to_creator[schema_idx]->Create(is_variable_length_field);
    }
  private:
    std::vector<std::shared_ptr<VariantFieldCreatorBase>> m_schema_idx_to_creator;
};
#endif
