#ifndef VARIANT_CELL_H
#define VARIANT_CELL_H

#include "headers.h"
#include "array_schema.h"

class VariantQueryConfig;
/*
 * This class is useful for storing cell data where all fields are 
 * stored in a single contiguous buffer
 */
class BufferVariantCell
{
  public:
    //Iterator over fields
    class FieldsIter
    {
      public:
       FieldsIter(const BufferVariantCell* ptr, const size_t idx)
       {
         m_ptr = ptr;
         m_idx = idx;
       }
       inline bool operator!=(const FieldsIter& other) { return (other.m_ptr != m_ptr || other.m_idx != m_idx); }
       template<typename T=void>
       inline const T* operator*() const { return m_ptr->get_field_ptr_for_query_idx<T>(m_idx); }
       inline int get_field_length() const { return m_ptr->get_field_length(m_idx); }
       inline bool is_variable_length_field() const { return m_ptr->is_variable_length_field(m_idx); }
       inline const FieldsIter& operator++()
       {
         ++m_idx;
         return *this;
       }
      private:
        const BufferVariantCell* m_ptr;
        size_t m_idx;
    };
  public:
    BufferVariantCell(const ArraySchema& array_schema, const VariantQueryConfig& query_config);
    void clear();
    void set_cell(const void* ptr);
    template<typename T=void>
    inline const T* get_field_ptr_for_query_idx(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_field_ptrs.size());
      return reinterpret_cast<const T*>(m_field_ptrs[query_idx]);
    }
    inline int get_field_length(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_field_lengths.size());
      return m_field_lengths[query_idx];
    }
    inline int get_field_length_descriptor(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_field_lengths.size());
      return m_field_length_descriptors[query_idx];
    }
    inline bool is_variable_length_field(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_field_lengths.size());
      return (m_field_length_descriptors[query_idx] == VAR_SIZE);
    }
    FieldsIter begin() const { return FieldsIter(this, 0ull); }
    FieldsIter end() const { return FieldsIter(this, m_field_ptrs.size()); }
  private:
    const ArraySchema* m_array_schema;
    const VariantQueryConfig* m_query_config;
    const uint8_t* m_cell_ptr;
    std::vector<const void*> m_field_ptrs;
    std::vector<size_t> m_field_element_sizes;
    std::vector<int> m_field_length_descriptors;
    std::vector<int> m_field_lengths;
};

#endif
