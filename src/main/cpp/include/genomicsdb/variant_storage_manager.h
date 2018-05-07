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

#ifndef VARIANT_STORAGE_MANAGER_H
#define VARIANT_STORAGE_MANAGER_H

#include "headers.h"
#include "variant_array_schema.h"
#include "variant_cell.h"
#include "genomicsdb_iterators.h"
#include <zlib.h>
#include "tiledb.h"
#include "timer.h"

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
#ifdef DO_PROFILING
      m_tiledb_timer.print("TileDB iterator", std::cerr);
      m_tiledb_to_buffer_cell_timer.print("TileDB to buffer cell", std::cerr);
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
#ifdef DO_PROFILING
      m_tiledb_timer.start();
#endif
      auto status = tiledb_array_iterator_next(m_tiledb_array_iterator);
      if(status != TILEDB_OK)
        throw VariantStorageManagerException(std::string("VariantArrayCellIterator increment failed")
            + "\nTileDB error message : "+tiledb_errmsg);
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
#ifdef DO_PROFILING
      m_tiledb_timer.stop();
#endif
      return *this;
    }
    const BufferVariantCell& operator*();
    //Set new interval to query
    void reset_subarray(const int64_t* coords)
    {
      assert(m_tiledb_array_iterator);
      if(tiledb_array_iterator_reset_subarray(m_tiledb_array_iterator,
          reinterpret_cast<const void*>(coords)) != TILEDB_OK)
        throw VariantStorageManagerException("Error resetting TileDB iterator subarray");
    }
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
#ifdef DO_PROFILING
    Timer m_tiledb_timer;
    Timer m_tiledb_to_buffer_cell_timer;
#endif
};

class VariantArrayInfo
{
  public:
    VariantArrayInfo(int idx, int mode,
        const std::string& workspace, const std::string& name,
        const VidMapper* vid_mapper,
        const VariantArraySchema& schema,
        TileDB_CTX* tiledb_ctx,
        TileDB_Array* tiledb_array,
        const size_t buffer_size=10u*1024u*1024u); //10MB buffer
    //Delete default copy constructor as it is incorrect
    VariantArrayInfo(const VariantArrayInfo& other) = delete;
    //Define move constructor explicitly
    VariantArrayInfo(VariantArrayInfo&& other);
    ~VariantArrayInfo()
    {
      close_array();
    }
    void close_array(const bool consolidate_tiledb_array=false);
    void set_schema(const VariantArraySchema& schema)
    {
      m_schema = schema;
      m_cell = std::move(BufferVariantCell(m_schema));
    }
    const VariantArraySchema& get_schema() const { return m_schema; }
    const std::string& get_array_name() const { return m_name; }
    void write_cell(const void* ptr);
    //Read #valid rows from metadata if available, else set from schema (array domain)
    void read_row_bounds_from_metadata();
    int read_row_bounds_from_metadata(const std::string& filename, const std::string& filepath,
        const unsigned char filetype);
    /*
     * Update #valid rows in the metadata
     */
    void update_row_bounds_in_array(
        const int64_t lb_row_idx, const int64_t max_valid_row_idx_in_array);
    //Return #valid rows in the array
    inline int64_t get_num_valid_rows_in_array() const
    {
      return (m_max_valid_row_idx_in_array - m_schema.dim_domains()[0].first + 1);
    }
    inline const VidMapper* get_vid_mapper() const { return m_vid_mapper; }
    const TileDB_Array* get_tiledb_array() const { return m_tiledb_array; }
  private:
    void get_genomicsdb_meta_directory_entries(std::string& buffer,
        std::vector<char*>& dir_entries_pointers,
        std::vector<unsigned char>& dir_entries_type);
  private:
    int m_idx;
    int m_mode;
    std::string m_workspace;
    std::string m_name;
    VariantArraySchema m_schema;
    BufferVariantCell m_cell;
    TileDB_CTX* m_tiledb_ctx;
    TileDB_Array* m_tiledb_array;
    //For writing cells
    //Buffers to hold data
    std::vector<std::vector<uint8_t>> m_buffers;
    //Pointers to buffers
    std::vector<void*> m_buffer_pointers;
    //Buffer offsets - byte where next data item needs to be written
    std::vector<size_t> m_buffer_offsets;
    //Max valid row idx in array
    int64_t m_max_valid_row_idx_in_array;
    bool m_metadata_contains_max_valid_row_idx_in_array;
    //Lb row idx
    bool m_metadata_contains_lb_row_idx;
    int64_t m_lb_row_idx;
    //Vid mapper
    const VidMapper* m_vid_mapper;
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
    }
    //Delete move and copy constructors
    VariantStorageManager(const VariantStorageManager& other) = delete;
    VariantStorageManager(VariantStorageManager&& other) = delete;
    /*
     * Wrapper functions around the C-API
     */
    bool check_if_TileDB_array_exists(const std::string& array_name);
    int open_array(const std::string& array_name, const VidMapper* vid_mapper, const char* mode,
        const int tiledb_zlib_compression_level=TILEDB_COMPRESSION_LEVEL_GZIP);
    void close_array(const int ad, const bool consolidate_tiledb_array=false);
    int define_array(const VariantArraySchema* variant_array_schema, const size_t num_cells_per_tile=1000u);
    void delete_array(const std::string& array_name);
    int define_metadata_schema(const VariantArraySchema* variant_array_schema);
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
     * For columnar iterator
     */
    SingleCellTileDBIterator* begin_columnar_iterator(
        int ad, const VariantQueryConfig& query_config,
        const bool use_common_array_object) const;
    /*
     * Write sorted cell
     */
    void write_cell_sorted(const int ad, const void* ptr);
    /*
     * Return #valid rows in the array
     */
    int64_t get_num_valid_rows_in_array(const int ad) const;
    /*
     * Update row bounds in the metadata
     */
    void update_row_bounds_in_array(const int ad, const int64_t lb_row_idx, const int64_t max_valid_row_idx_in_array);
    /*
     * Return workspace path
     */
    const std::string& get_workspace() const { return m_workspace; }
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
    //Metadata attribute name
    static std::vector<const char*> m_metadata_attributes;
};

#endif
