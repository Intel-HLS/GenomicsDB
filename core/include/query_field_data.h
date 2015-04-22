#ifndef QUERY_FIELD_DATA
#define QUERY_FIELD_DATA

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <memory>
#include "array_schema.h"
#include <typeindex>
#include <unordered_map>

class UnknownAttributeTypeException {
  public:
    UnknownAttributeTypeException(const std::string m="Unknown type of queried attribute") : msg_(m) { ; }
    ~UnknownAttributeTypeException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const std::string& what() const { return msg_; }
  private:
    std::string msg_;
};

class QueryFieldData
{
  public:
    /**
     * Do nothing constructor
     */
    QueryFieldData();
    /**
     * Destructor that just verifies that data ptr was not allocated
     */
    ~QueryFieldData()
    {
      assert(!m_is_data_ptr_allocated); //should have been explicitly de-allocated, if allocated earlier
    }
    void clear() { ; }
    void reset();
    /**
     * Return an element of the data with correct casting
     */
    template<class T>
    inline T get_element(uint64_t idx)
    {
      assert(idx < m_num_elements);
      return (static_cast<T*>(m_data_ptr))[idx];
    }
    /**
     * Return pointer to the data with the right casting at element idx
     */
    template<class T>
    inline T* get_pointer(uint64_t idx)
    {
      return (static_cast<T*>(m_data_ptr))+idx;
    }
    /**
     * Return pointer to the data with the right casting at element 0
     */
    template<class T>
    T* get_pointer() { return get_pointer<T>(0ull); }
    /**
     * Set field info
     */
    inline void set_field_info(const std::type_index& element_type)
    {
      m_element_type = element_type;
      const auto iter = QueryFieldData::m_type_to_size.find(element_type);
      if(iter == QueryFieldData::m_type_to_size.end())
        throw UnknownAttributeTypeException("Unknown type of queried attribute "+std::string(element_type.name()));
      else
        m_element_size = (*iter).second;
    }
    /**
     * Allocate memory using Allocator
     */
    template<class AllocatorTy>
    inline void allocate_memory(uint64_t num_elements, AllocatorTy A)
    {
      m_data_ptr = static_cast<void*>(A.allocate(num_elements));
      m_num_elements = num_elements;
      m_is_data_ptr_allocated = true;
    }
    /**
     * De-allocate memory using Allocator
     */
    template<class AllocatorTy>
    inline void deallocate_memory(AllocatorTy A)
    {
      assert(m_is_data_ptr_allocated);
      typedef typename std::allocator_traits<AllocatorTy>::pointer PtrTy;       //ptr type of real data
      A.deallocate(static_cast<PtrTy>(m_data_ptr), m_num_elements);
      m_is_data_ptr_allocated = false;
    }
    /**
     * Copy data internally, memory MUST be allocated before calling this function
     */
    inline void copy_data(void* src, uint64_t num_elements, uint64_t offset_in_bytes=0ull)
    {
      assert(m_is_data_ptr_allocated && offset_in_bytes+num_elements*m_element_size <= m_num_elements*m_element_size);
      memcpy(static_cast<char*>(m_data_ptr)+offset_in_bytes, src, num_elements*m_element_size);
    }
    /**
     * Don't copy, but make data_ptr point to this address
     */
    inline void point_to(void* src, uint64_t num_elements)
    {
      assert(!m_is_data_ptr_allocated);
      m_data_ptr = src;
      m_num_elements = num_elements;
    }
  private:
    void* m_data_ptr;           //ptr to data
    /** The following variable is:
     * true if the data_ptr points to a memory segment allocated by this object
     * false if the ptr simply points to a region of memory which is managed elsewhere
     */
    bool m_is_data_ptr_allocated;
    std::type_index m_element_type;
    size_t m_element_size;    //sizeof(float/int etc)
    uint64_t m_num_elements;    //all data is a vector - single elements = vector of size 1
    /*
     * Static member that maps std::type_index to size
     */
    static std::unordered_map<std::type_index, size_t> m_type_to_size;
};

#endif
