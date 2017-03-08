/**
 * The MIT License (MIT)
 * Copyright (c) 2016 Intel Corporation
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

#ifndef VARIANT_CELL_H
#define VARIANT_CELL_H

#include "headers.h"
#include "variant_array_schema.h"

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
    BufferVariantCell(const VariantArraySchema& array_schema);
    BufferVariantCell(const VariantArraySchema& array_schema, const VariantQueryConfig& query_config);
    BufferVariantCell(const VariantArraySchema& array_schema, const std::vector<int>& attribute_ids);
    void clear();
    void set_variant_array_schema(const VariantArraySchema& array_schema)
    {
      m_array_schema = &array_schema;
    }
    void set_cell(const void* ptr);
    template<typename T=void>
    inline const T* get_field_ptr_for_query_idx(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_field_ptrs.size());
      return reinterpret_cast<const T*>(m_field_ptrs[query_idx]);
    }
    inline void set_field_ptr_for_query_idx(const int query_idx, const void* ptr)
    {
      assert(static_cast<size_t>(query_idx) < m_field_ptrs.size());
      m_field_ptrs[query_idx] = ptr;
    }
    inline int get_field_length(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_field_lengths.size());
      return m_field_lengths[query_idx];
    }
    inline size_t get_field_size_in_bytes(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_schema_idxs.size());
      auto schema_idx = m_schema_idxs[query_idx];
      return m_field_lengths[query_idx]*(m_array_schema->element_size(schema_idx));
    }
    inline void set_field_length(const int query_idx, const int length)
    {
      assert(static_cast<size_t>(query_idx) < m_field_lengths.size());
      m_field_lengths[query_idx] = length;
    }
    inline void set_field_size_in_bytes(const int query_idx, const size_t bytes)
    {
      assert(static_cast<size_t>(query_idx) < m_schema_idxs.size());
      auto schema_idx = m_schema_idxs[query_idx];
      m_field_lengths[query_idx] = bytes/(m_array_schema->element_size(schema_idx));
    }
    inline int get_field_length_descriptor(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_schema_idxs.size());
      auto schema_idx = m_schema_idxs[query_idx];
      return m_array_schema->val_num(schema_idx);
    }
    inline bool is_variable_length_field(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_schema_idxs.size());
      auto schema_idx = m_schema_idxs[query_idx];
      return m_array_schema->is_variable_length_field(schema_idx);
    }
    FieldsIter begin() const { return FieldsIter(this, 0ull); }
    FieldsIter end() const { return FieldsIter(this, m_field_ptrs.size()); }
    inline int64_t get_begin_column() const { return m_begin_column_idx; }
    inline int64_t get_row() const { return m_row_idx; }
    inline void set_coordinates(const int64_t row, const int64_t begin_column)
    {
      m_row_idx = row;
      m_begin_column_idx = begin_column;
    }
  private:
    void resize(const size_t num_fields);
    void update_field_info(const int query_idx, const int schema_idx);
  private:
    const VariantArraySchema* m_array_schema;
    std::vector<int> m_schema_idxs;
    std::vector<const void*> m_field_ptrs;
    std::vector<int> m_field_lengths;
    //Co-ordinates
    int64_t m_row_idx;
    int64_t m_begin_column_idx;
};

#endif
