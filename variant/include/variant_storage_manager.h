#ifndef VARIANT_STORAGE_MANAGER_H
#define VARIANT_STORAGE_MANAGER_H

#include "headers.h"
#include "storage_manager.h"
#include "variant_array_schema.h"
/*
 * Wrapper class around StorageManager - shields GenomicsDB from changes in
 * core
 */
class VariantStorageManager
{
  public:
    VariantStorageManager(const std::string& workspace, const unsigned segment_size=10u*1024u*1024u)
      : m_storage_manager(workspace, segment_size)
    {}
    /*
     * Wrapper functions around StorageManager
     */
    int open_array(const std::string& array_name, const char* mode)
    {
      return m_storage_manager.open_array(array_name, mode);
    }
    void close_array(const int ad) { m_storage_manager.close_array(ad); }
    int define_array(const VariantArraySchema* variant_array_schema)
    {
      return m_storage_manager.define_array(variant_array_schema->get_array_schema());
    }
    int get_array_schema(const int ad, VariantArraySchema* variant_array_schema) const
    {
      const ArraySchema* tmp_schema = 0;
      auto status = m_storage_manager.get_array_schema(ad, tmp_schema);
      if(status == TILEDB_OK)
        *variant_array_schema = std::move(VariantArraySchema(tmp_schema));
      return status;
    }
    /*
     * Wrapper around forward iterator
     */
    ArrayConstCellIterator<int64_t>* begin(
        int ad, const int64_t* range, const std::vector<int>& attribute_ids) const 
    {
      return m_storage_manager.begin<int64_t>(ad, range, attribute_ids);
    }
    /*
     * Wrapper around reverse iterator
     */
    ArrayConstReverseCellIterator<int64_t>* rbegin(
        int ad, const int64_t* range, const std::vector<int>& attribute_ids,
        bool is_range = true) const 
    {
      return m_storage_manager.rbegin<int64_t>(ad, range, attribute_ids, is_range);
    }
    /*
     * Write sorted cell
     */
    void write_cell_sorted(const int ad, const void* cell_ptr)
    {
      m_storage_manager.write_cell_sorted<int64_t>(ad, cell_ptr);
    }
  private:
    StorageManager m_storage_manager;
};

#endif
