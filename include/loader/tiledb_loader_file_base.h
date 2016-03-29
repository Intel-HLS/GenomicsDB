#ifndef TILEDB_IMPORTER_H
#define TILEDB_IMPORTER_H

#include "gt_common.h"

class File2TileDBBinaryColumnPartitionBase
{
  public:
    File2TileDBBinaryColumnPartitionBase() { ; }
  protected:
    int64_t m_column_interval_begin;
    int64_t m_column_interval_end;
    //Buffer offsets - 1 per callset
    //Offset at which data should be copied for the current batch
    std::vector<int64_t> m_begin_buffer_offset_for_local_callset;
    //Offset at which the current line begins
    std::vector<int64_t> m_last_full_line_end_buffer_offset_for_local_callset;
    //Current value of offset
    std::vector<int64_t> m_buffer_offset_for_local_callset;
};

/*
 * Base class for all instances of objects that convert file formats for importing
 * to TileDB
 */

class File2TileDBBinaryBase
{
  public:
    File2TileDBBinaryBase(const std::string& filename,
        unsigned file_idx, VidMapper& vid_mapper,
        size_t max_size_per_callset,
        bool treat_deletions_as_intervals,
        bool parallel_partitions=false, bool noupdates=true, bool close_file=false)
    {
      m_filename = filename;
      m_file_idx = file_idx;
      m_vid_mapper = &(vid_mapper);
      m_max_size_per_callset = max_size_per_callset;
      m_treat_deletions_as_intervals = treat_deletions_as_intervals;
      m_parallel_partitions = parallel_partitions;
      m_noupdates = noupdates;
      m_close_file = close_file;
    }
    void clear()
    {
      m_filename.clear();
      m_local_callset_idx_to_tiledb_row_idx.clear();
      m_enabled_local_callset_idx_vec.clear();
    }
  protected:
    bool m_parallel_partitions;
    bool m_noupdates;
    bool m_close_file;
    bool m_treat_deletions_as_intervals;
    VidMapper* m_vid_mapper;
    int64_t m_file_idx;
    size_t m_max_size_per_callset;
    std::string m_filename;
    //Local callset idx to tiledb row idx
    std::vector<int64_t> m_local_callset_idx_to_tiledb_row_idx;
    //Enabled local callset idx
    std::vector<int64_t> m_enabled_local_callset_idx_vec;
};

#endif
