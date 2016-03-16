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
    }
    inline bool end() const {
      return tiledb_array_iterator_end(m_tiledb_array_iterator);
    }
    inline const VariantArrayCellIterator& operator++()
    {
      tiledb_array_iterator_next(m_tiledb_array_iterator);
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
};

class ArrayInfo
{
  public:
    int m_idx;
    int m_mode;
    std::string m_name;
    const VariantArraySchema* m_schema;
};

/*
 * Wrapper class around TileDB C API - shields GenomicsDB from changes in
 * core
 */
class VariantStorageManager
{
  public:
    VariantStorageManager(const std::string& workspace, const unsigned segment_size=10u*1024u*1024u)
    {
      m_workspace = workspace;
      m_segment_size = segment_size;
      /*Initialize context with default params*/
      tiledb_ctx_init(&m_tiledb_ctx, NULL);
    }
    ~VariantStorageManager()
    {
      m_open_arrays_info_vector.clear();
      m_workspace.clear();
       /* Finalize context. */
      tiledb_ctx_finalize(m_tiledb_ctx);
    }
    /*
     * Wrapper functions around the C-API
     */
    int open_array(const std::string& array_name, const char* mode);
    void close_array(const int ad);
    int define_array(const VariantArraySchema* variant_array_schema);
    int get_array_schema(const int ad, VariantArraySchema* variant_array_schema);
    /*
     * Wrapper around forward iterator
     */
    VariantArrayCellIterator begin(
        int ad, const int64_t* range, const std::vector<int>& attribute_ids) const ;
    /*
     * Write sorted cell
     * FIXME
     */
    void write_cell_sorted(const int ad, const void* cell_ptr) { ; }
  private:
    static const std::unordered_map<std::string, int> m_mode_string_to_int;
    //TileDB context
    TileDB_CTX* m_tiledb_ctx;
    //Workspace name
    std::string m_workspace;
    //Info vector for open arrays
    std::vector<ArrayInfo> m_open_arrays_info_vector;
    //How much data to read in a given access
    size_t m_segment_size;
};

#endif
