/**
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

#ifndef VARIANT_STORAGE_MANAGER_H
#define VARIANT_STORAGE_MANAGER_H

#include "headers.h"
#include "variant_array_schema.h"
#include "variant_cell.h"
#include "c_api.h"

//Exceptions thrown 
class VariantStorageManagerException : public std::exception {
  public:
    VariantStorageManagerException(const std::string m="") : msg_("VariantStorageManagerException exception : "+m) { ; }
    ~VariantStorageManagerException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class VariantArrayCellIterator
{
  public:
    VariantArrayCellIterator(TileDB_CTX* tiledb_ctx, const VariantArraySchema& variant_array_schema,
        const std::string& array_path, const int64_t* range, const std::vector<int>& attribute_ids, const size_t buffer_size);
    ~VariantArrayCellIterator()
    {
      if(m_tiledb_array_iterator)
        tiledb_array_iterator_finalize(m_tiledb_array_iterator);
      m_tiledb_array_iterator = 0;
#ifdef DEBUG
      std::cerr << "#cells traversed "<<m_num_cells_iterated_over<<"\n";
#endif
    }
    //Delete copy and move constructors
    VariantArrayCellIterator(const VariantArrayCellIterator& other) = delete;
    VariantArrayCellIterator(VariantArrayCellIterator&& other) = delete;
    inline bool end() const {
      return tiledb_array_iterator_end(m_tiledb_array_iterator);
    }
    inline const VariantArrayCellIterator& operator++()
    {
      tiledb_array_iterator_next(m_tiledb_array_iterator);
#ifdef DEBUG
      if(!end())
      {
        ++m_num_cells_iterated_over;
        //Co-ordinates
        const uint8_t* field_ptr = 0;
        size_t field_size = 0u;
        tiledb_array_iterator_get_value(m_tiledb_array_iterator, m_num_queried_attributes,
            reinterpret_cast<const void**>(&field_ptr), &field_size);
        assert(field_size == m_variant_array_schema->dim_size_in_bytes());
        auto coords_ptr = reinterpret_cast<const int64_t*>(field_ptr);
        assert(coords_ptr[1] > m_last_column || (coords_ptr[1] == m_last_column && coords_ptr[0] > m_last_row));
        m_last_row = coords_ptr[0];
        m_last_column = coords_ptr[1];
      }
#endif
      return *this;
    }
    const BufferVariantCell& operator*();
  private:
    unsigned m_num_queried_attributes;
    TileDB_CTX* m_tiledb_ctx;
    const VariantArraySchema* m_variant_array_schema;
    BufferVariantCell m_cell;
    //The actual TileDB array iterator
    TileDB_ArrayIterator* m_tiledb_array_iterator;
    //Buffers to hold data
    std::vector<std::vector<uint8_t>> m_buffers;
    //Pointers to buffers
    std::vector<const void*> m_buffer_pointers;
    //Buffer sizes
    std::vector<size_t> m_buffer_sizes;
#ifdef DEBUG
    int64_t m_last_row;
    int64_t m_last_column;
    uint64_t m_num_cells_iterated_over;
#endif
};

class VariantArrayInfo
{
  public:
    VariantArrayInfo(int idx, int mode, const std::string& name, const VariantArraySchema& schema, TileDB_Array* tiledb_array,
        const size_t buffer_size=10u*1024u*1024u); //10MB buffer
    //Delete default copy constructor as it is incorrect
    VariantArrayInfo(const VariantArrayInfo& other) = delete;
    //Define move constructor explicitly
    VariantArrayInfo(VariantArrayInfo&& other);
    ~VariantArrayInfo()
    {
      close_array();
    }
    void close_array()
    {
      //Flush cells in buffer
      auto coords_buffer_idx = m_buffers.size()-1u;
      if((m_mode == TILEDB_ARRAY_WRITE || m_mode == TILEDB_ARRAY_WRITE_UNSORTED)
          && m_buffer_offsets[coords_buffer_idx] > 0ull)
      {
        tiledb_array_write(m_tiledb_array, const_cast<const void**>(&(m_buffer_pointers[0])), &(m_buffer_offsets[0]));
        memset(&(m_buffer_offsets[0]), 0, m_buffer_offsets.size()*sizeof(size_t));
      }
      if(m_tiledb_array)
        tiledb_array_finalize(m_tiledb_array);
      m_tiledb_array = 0;
      m_name.clear();
      m_mode = -1;
    }
    void set_schema(const VariantArraySchema& schema)
    {
      m_schema = schema;
      m_cell = std::move(BufferVariantCell(m_schema));
    }
    const VariantArraySchema& get_schema() const { return m_schema; }
    const std::string& get_array_name() const { return m_name; }
    void write_cell(const void* ptr);
  private:
    int m_idx;
    int m_mode;
    std::string m_name;
    VariantArraySchema m_schema;
    BufferVariantCell m_cell;
    TileDB_Array* m_tiledb_array;
    //For writing cells
    //Buffers to hold data
    std::vector<std::vector<uint8_t>> m_buffers;
    //Pointers to buffers
    std::vector<void*> m_buffer_pointers;
    //Buffer offsets - byte where next data item needs to be written
    std::vector<size_t> m_buffer_offsets;
#ifdef DEBUG
    int64_t m_last_row;
    int64_t m_last_column;
#endif
};

/*
 * Wrapper class around TileDB C API - shields GenomicsDB from changes in
 * core
 */
class VariantStorageManager
{
  public:
    VariantStorageManager(const std::string& workspace, const unsigned segment_size=10u*1024u*1024u);
    ~VariantStorageManager()
    {
      m_open_arrays_info_vector.clear();
      m_workspace.clear();
       /* Finalize context. */
      tiledb_ctx_finalize(m_tiledb_ctx);
      free(m_tiledb_ctx);
    }
    //Delete move and copy constructors
    VariantStorageManager(const VariantStorageManager& other) = delete;
    VariantStorageManager(VariantStorageManager&& other) = delete;
    /*
     * Wrapper functions around the C-API
     */
    int open_array(const std::string& array_name, const char* mode);
    void close_array(const int ad);
    int define_array(const VariantArraySchema* variant_array_schema, const size_t num_cells_per_tile=1000u);
    /*
     * Load array schema
     */
    int get_array_schema(const std::string& array_name, VariantArraySchema* variant_array_schema);
    int get_array_schema(const int ad, VariantArraySchema* variant_array_schema);
    /*
     * Wrapper around forward iterator
     */
    VariantArrayCellIterator* begin(
        int ad, const int64_t* range, const std::vector<int>& attribute_ids) const ;
    /*
     * Write sorted cell
     */
    void write_cell_sorted(const int ad, const void* ptr);
  private:
    static const std::unordered_map<std::string, int> m_mode_string_to_int;
    //TileDB context
    TileDB_CTX* m_tiledb_ctx;
    //Workspace name
    std::string m_workspace;
    //Info vector for open arrays
    std::vector<VariantArrayInfo> m_open_arrays_info_vector;
    //How much data to read/write in a given access
    size_t m_segment_size;
};

#endif
