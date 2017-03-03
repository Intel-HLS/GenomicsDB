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
#include "variant_cell.h"
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

/*
 * Iterates over TileDB cells one at a time
 */
class SingleCellTileDBIterator
{
  public:
    SingleCellTileDBIterator(TileDB_CTX* tiledb_ctx, const VariantArraySchema& variant_array_schema,
        const std::string& array_path, const int64_t* range, const std::vector<int>& attribute_ids, const size_t buffer_size);
    ~VariantArrayCellIterator()
    {
      if(m_tiledb_array_iterator)
        tiledb_array_iterator_finalize(m_tiledb_array_iterator);
      m_tiledb_array_iterator = 0;
#ifdef DO_PROFILING
      m_tiledb_timer.print("TileDB iterator", std::cerr);
      m_tiledb_to_buffer_cell_timer.print("TileDB to buffer cell", std::cerr);
#endif
    }
    //Delete copy and move constructors
    VariantArrayCellIterator(const VariantArrayCellIterator& other) = delete;
    VariantArrayCellIterator(VariantArrayCellIterator&& other) = delete;
    //Iterator functionality
    const BufferVariantCell& operator*();
    inline const VariantArrayCellIterator& operator++();
    inline bool end() const {
      return tiledb_array_iterator_end(m_tiledb_array_iterator);
    }
  private:
    TileDB_CTX* m_tiledb_ctx;
    const VariantArraySchema* m_variant_array_schema;
    BufferVariantCell m_cell;
    //Buffers for fields
    std::vector<GenomicsDBColumnarField> m_fields;
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
