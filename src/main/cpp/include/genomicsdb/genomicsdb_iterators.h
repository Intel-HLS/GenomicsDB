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

#ifndef GENOMICSDB_ITERATORS_H
#define GENOMICSDB_ITERATORS_H

#include "headers.h"
#include "variant_array_schema.h"
#include "genomicsdb_columnar_field.h"
#include "c_api.h"
#include "timer.h"

//Exceptions thrown
class GenomicsDBIteratorException : public std::exception {
  public:
    GenomicsDBIteratorException(const std::string m="") : msg_("GenomicsDBIteratorException exception : "+m) { ; }
    ~GenomicsDBIteratorException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class GenomicsDBColumnarCell;
/*
 * Iterates over TileDB cells one at a time
 */
class SingleCellTileDBIterator
{
  public:
    SingleCellTileDBIterator(TileDB_CTX* tiledb_ctx, const VariantArraySchema& variant_array_schema,
        const std::string& array_path, const int64_t* range, const std::vector<int>& attribute_ids, const size_t buffer_size);
    ~SingleCellTileDBIterator();
    //Delete copy and move constructors
    SingleCellTileDBIterator(const SingleCellTileDBIterator& other) = delete;
    SingleCellTileDBIterator(SingleCellTileDBIterator&& other) = delete;
    //Iterator functionality
    inline const GenomicsDBColumnarCell& operator*() const
    {
      return *m_cell;
    }
    inline const SingleCellTileDBIterator& operator++();
    //Get field pointer, length, size
    inline const uint8_t* get_field_ptr_for_query_idx(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_fields.size());
      auto& genomicsdb_columnar_field = m_fields[query_idx];
      return genomicsdb_columnar_field.get_pointer_to_curr_index_data_in_live_list_tail();
    }
    inline int get_field_length(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_fields.size());
      auto& genomicsdb_columnar_field = m_fields[query_idx];
      return genomicsdb_columnar_field.get_length_of_curr_index_data_in_live_list_tail();
    }
    inline size_t get_field_size_in_bytes(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_fields.size());
      auto& genomicsdb_columnar_field = m_fields[query_idx];
      return genomicsdb_columnar_field.get_size_of_curr_index_data_in_live_list_tail();
    }
    inline bool is_valid(const int query_idx) const
    {
      assert(static_cast<size_t>(query_idx) < m_fields.size());
      auto& genomicsdb_columnar_field = m_fields[query_idx];
      return genomicsdb_columnar_field.is_valid_curr_index_data_in_live_list_tail();
    }
  protected:
    /*
     * Does one read for the attributes in m_query_attribute_idx_vec
     */
    void read_from_TileDB();
  private:
    const VariantArraySchema* m_variant_array_schema;
    GenomicsDBColumnarCell* m_cell;
    //Buffers for fields
    std::vector<GenomicsDBColumnarField> m_fields;
    //Contains query idx for only the fields that must be fetched from TileDB in the next round
    //The first time all fields are queried - in subsequent iterations only those fields whose
    //buffers are consumed completely are queried
    std::vector<int> m_query_attribute_idx_vec;
    //Since variable length fields have buffers for offsets, need a mapping structure
    std::vector<size_t> m_query_attribute_idx_to_tiledb_buffer_idx;
    std::vector<void*> m_buffer_pointers;
    std::vector<size_t> m_buffer_sizes;
    //The TileDB array object
    TileDB_Array* m_tiledb_array;
#ifdef DO_PROFILING
    Timer m_tiledb_timer;
    Timer m_tiledb_to_buffer_cell_timer;
#endif
};

#endif
