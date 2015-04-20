#ifndef QUERY_FIELD_DATA
#define QUERY_FIELD_DATA

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <memory>
#include "array_schema.h"

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
    inline void set_field_info(unsigned schema_idx, const std::type_info* element_type)
    {
      m_schema_idx = schema_idx;
      m_element_type = element_type;
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
    unsigned m_schema_idx;      //attribute idx in the array
    const std::type_info* m_element_type;
    unsigned m_element_size;    //sizeof(float/int etc)
    uint64_t m_num_elements;    //all data is a vector - single elements = vector of size 1
};

#endif
