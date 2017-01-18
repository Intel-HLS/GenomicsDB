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

#ifndef GENOMICSDB_COLUMNAR_FIELD_H
#define GENOMICSDB_COLUMNAR_FIELD_H

#include "headers.h"
#include "vcf.h"

class GenomicsDBColumnarFieldException : public std::exception {
  public:
    GenomicsDBColumnarFieldException(const std::string m) : msg_(m) { ; }
    ~GenomicsDBColumnarFieldException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

/*
 * Think of this class as a vector which cannot be resized
 * The buffer space of this object is owned by GenomicsDBColumnarField - this object
 * only contains offsets within the buffer
 */
template<class DataType>
class LinkedVector
{

};

/*
 * Class that stores data for a given field from multiple VariantCalls/TileDB cells together in a columnar fashion
 * Maintains a big buffer for data, valid bitvector and offsets array
 */
class GenomicsDBColumnarField
{
  public:
    GenomicsDBColumnarField(const std::type_index element_type, const int length_descriptor,
        const unsigned fixed_length_field_length);
    void clear()
    {
      m_data.clear();
      m_valid.clear();
      m_begin_offsets.clear();
    }
    void resize(const uint64_t num_cells);
    /*
     * Put data into the columnar structure
     * Check validity of structure
     */
    void put(const uint64_t cell_idx, const uint8_t* src, const size_t num_bytes)
    {
      assert(cell_idx < m_valid.size());
      auto is_valid = m_check_validity(src, (m_length_descriptor == BCF_VL_FIXED) ? m_fixed_length_field_size : num_bytes);
      m_valid[cell_idx] = is_valid;
      if(is_valid)
        if(m_length_descriptor == BCF_VL_FIXED)
        {
          assert(m_fixed_length_field_size*cell_idx+m_fixed_length_field_size <= m_data.size());
          memcpy(&(m_data[m_fixed_length_field_size*cell_idx]), src, m_fixed_length_field_size);
        }
        else
        {
          resize(cell_idx, num_bytes);
          memcpy(&(m_data[m_begin_offsets[cell_idx]]), src, num_bytes);
          m_num_bytes[cell_idx] = num_bytes;
        }
      else
        if(m_length_descriptor != BCF_VL_FIXED)
        {
          assert(cell_idx < m_num_bytes.size());
          m_num_bytes[cell_idx] = 0u;
        }
    }
    void resize(const uint64_t cell_idx, const size_t num_bytes)
    {
      if(m_length_descriptor == BCF_VL_FIXED)
        return;
      assert(cell_idx+1u < m_begin_offsets.size()); //see sizing of m_begin_offsets
      assert(m_begin_offsets[cell_idx+1] >= m_begin_offsets[cell_idx]+m_num_bytes[cell_idx]);
      if(m_begin_offsets[cell_idx]+num_bytes > m_begin_offsets[cell_idx+1])
      {
        auto extra = (m_begin_offsets[cell_idx]+num_bytes-m_begin_offsets[cell_idx+1]);
        auto old_size = m_data.size();
        m_data.resize(m_data.size() + extra);
        memmove(&(m_data[m_begin_offsets[cell_idx]+num_bytes]), &(m_data[m_begin_offsets[cell_idx+1]]),
            old_size-m_begin_offsets[cell_idx+1]);
        for(auto i=cell_idx+1u;i<m_begin_offsets.size();++i)
          m_begin_offsets[i] += extra;
      }
    }
    /*
     * Is valid?
     */
    bool is_valid(const uint64_t cell_idx) const
    {
      assert(cell_idx < m_valid.size());
      return m_valid[cell_idx];
    }
    /*
     * Get pointer to data and number of bytes - must be valid (caller must check before calling)
     * Insertions may invalidate the pointer returned because of resizing of the data buffer
     */
    uint8_t* get(const uint64_t cell_idx)
    {
      assert(is_valid(cell_idx));
      if(m_length_descriptor == BCF_VL_FIXED)
        return &(m_data[m_fixed_length_field_size*cell_idx]);
      else
      {
        assert(cell_idx+1u < m_begin_offsets.size());
        return &(m_data[m_begin_offsets[cell_idx]]);
      }
    }
    const uint8_t* get(const uint64_t cell_idx) const
    {
      assert(is_valid(cell_idx));
      if(m_length_descriptor == BCF_VL_FIXED)
        return &(m_data[m_fixed_length_field_size*cell_idx]);
      else
      {
        assert(cell_idx+1u < m_begin_offsets.size());
        return &(m_data[m_begin_offsets[cell_idx]]);
      }
    }
    size_t size(const uint64_t cell_idx) const
    {
      if(m_length_descriptor == BCF_VL_FIXED)
        return m_fixed_length_field_size;
      else
        return m_num_bytes[cell_idx];
    }
    size_t length(const uint64_t cell_idx) const
    {
      if(m_length_descriptor == BCF_VL_FIXED)
        return m_fixed_length_field_num_elements;
      else
        return m_num_bytes[cell_idx]/m_element_size;
    }
    template<typename T>
    static bool check_validity(const uint8_t* ptr, const size_t num_bytes)
    {
      auto data = reinterpret_cast<const T*>(ptr);
      auto num_elements = num_bytes/sizeof(T);
      for(auto i=0ull;i<num_elements;++i)
        if(is_bcf_valid_value<T>(data[i]))
          return true;
      return false;
    }
  private:
    void assign_function_pointers();
  private:
    int m_length_descriptor;
    unsigned m_fixed_length_field_num_elements;
    unsigned m_fixed_length_field_size;
    unsigned m_element_size;
    std::type_index m_element_type;
    std::vector<uint8_t> m_data;
    std::vector<bool> m_valid;
    std::vector<size_t> m_num_bytes;
    std::vector<size_t> m_begin_offsets;
    //Function pointer that determines validity check
    bool (*m_check_validity)(const uint8_t* data, const size_t num_bytes);
};

#endif
