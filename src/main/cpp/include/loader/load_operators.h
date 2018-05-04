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

#ifndef LOAD_OPERATORS_H
#define LOAD_OPERATORS_H

#include "vid_mapper.h"
#include "query_variants.h"
#include "broad_combined_gvcf.h" 
#include "variant_storage_manager.h"
#include "json_config.h"

struct CellPointersColumnMajorCompare
{
  bool operator()(const uint8_t* const a, const uint8_t* const b) const
  {
    auto* a_coords = reinterpret_cast<const int64_t* const>(a);
    auto* b_coords = reinterpret_cast<const int64_t* const>(b);
    return (a_coords[1] < b_coords[1] || (a_coords[1] == b_coords[1] && a_coords[0] < b_coords[0]));
  }
};

//Exceptions thrown
class LoadOperatorException : public std::exception {
  public:
    LoadOperatorException(const std::string m="") : msg_("LoadOperatorException : "+m) { ; }
    ~LoadOperatorException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class LoaderOperatorBase
{
  public:
    LoaderOperatorBase(
      const std::string& loader_config_file,
      const size_t num_callsets,
      const int partition_idx,
      const VidMapper* vid_mapper=0,
      const bool vid_mapper_file_required = true)
      : m_column_partition(0, INT64_MAX-1), m_row_partition(0, INT64_MAX-1)
    {
      m_crossed_column_partition_begin = false;
      m_first_cell = true;
      m_partition_idx = partition_idx;
      m_loader_json_config.set_vid_mapper_file_required(
        vid_mapper_file_required);
#ifdef DUPLICATE_CELL_AT_END
      m_cell_copies.resize(num_callsets, 0);
      m_last_end_position_for_row.resize(num_callsets, -1ll);
#endif
      VidMapper tmp_vid_mapper;
      if(vid_mapper)
        tmp_vid_mapper = *vid_mapper;
      //Parse loader JSON
      m_loader_json_config.read_from_file(loader_config_file, &tmp_vid_mapper, partition_idx);
      //Partitioning information
      m_row_partition = RowRange(0, m_loader_json_config.get_max_num_rows_in_array()-1);
      if(m_loader_json_config.is_partitioned_by_row())
      {
        m_row_partition = m_loader_json_config.get_row_partition(partition_idx);
        m_row_partition.second = std::min(m_row_partition.second, static_cast<int64_t>(m_loader_json_config.get_max_num_rows_in_array()-1));
      }
      else
        m_column_partition = m_loader_json_config.get_column_partition(partition_idx);
    }
    virtual ~LoaderOperatorBase() { ; }
    /*
     * Function that is called before the parallel section takes over, but inside the operational loop
     */
    virtual void pre_operate_sequential() { ; }
    /*
     * Virtual function that must be overridden by sub-classes
     */
    virtual void operate(const void* cell_ptr) = 0;
    /*
     * Handle intervals which begin before the column partition and span across.
     * This is needed to deal with the overlapping variants issue for VCFs (see wiki)
     * @return: true if the interval begins before the column partition
     */
    void handle_intervals_spanning_partition_begin(const int64_t row, const int64_t begin, const int64_t end,
        const size_t cell_size, const void* cell_ptr);
    /*
     * Throws exceptions if cell does not belong to current partition
     */
    void check_cell_coordinates(const int64_t row_idx, const int64_t column_begin, const int64_t column_end);
    /*
     * Returns true if output buffer for this operator is full and caller must "block"
     * till flush is completed
     */
    virtual bool overflow() const { return false; }
    /*
     * Called within parallel sections - useful if the output needs to be flushed by a thread
     * not in the critical path
     */
    virtual void flush_output() { ; }
    /*
     * Function that is called after the parallel section, but inside the operational loop
     */
    virtual void post_operate_sequential() { ; }
    /*
     * Called at the end - the argument is the column interval end limit
     */
    virtual void finish(const int64_t column_interval_end);
  protected:
    int m_partition_idx;
    ColumnRange m_column_partition;
    RowRange m_row_partition;
    bool m_crossed_column_partition_begin;
    bool m_first_cell;
#ifdef DUPLICATE_CELL_AT_END
    //Copy of cell buffers
    std::vector<uint8_t*> m_cell_copies;
    //End position of last cell seen for current row
    std::vector<int64_t> m_last_end_position_for_row;
#endif
    JSONLoaderConfig m_loader_json_config;
};

class LoaderArrayWriter : public LoaderOperatorBase
{
  public:
    LoaderArrayWriter(
      const VidMapper* id_mapper,
      const std::string& config_filename,
      int rank,
      const bool vid_mapper_file_required);
    virtual ~LoaderArrayWriter()
    {
      if(m_schema)
        delete m_schema;
      if(m_storage_manager)
        delete m_storage_manager;
#ifdef DUPLICATE_CELL_AT_END
      for(auto ptr : m_cell_copies)
        free(ptr);
      m_cell_copies.clear();
#endif
    }
    virtual void operate(const void* cell_ptr);
    virtual void finish(const int64_t column_interval_end);
  private:
    int m_array_descriptor;
    VariantArraySchema* m_schema;
    VariantStorageManager* m_storage_manager;
#ifdef DUPLICATE_CELL_AT_END
    /*
     * Function that writes top element from the PQ to disk
     * If the top element is a begin cell that spans multiple columns, create an END copy
     * and adds to PQ again
     */
    void write_top_element_to_disk();
    //Mimics behavior of program heap for cell copies - minimizes number of frees/reallocations
    std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> m_memory_manager;
    //For use in priority queue
    typedef struct
    {
      int64_t m_row;
      int64_t m_begin_column;
      int64_t m_end_column;
      size_t m_idx_in_cell_copies_vector;
    } CellWrapper;
    struct ColumnMajorCellCompareGT 
    {
      bool operator()(const CellWrapper& a, const CellWrapper& b) const
      {
        return ((a.m_begin_column > b.m_begin_column) || (a.m_begin_column == b.m_begin_column && a.m_row > b.m_row));
      }
    };
    //top() contains CellWrapper with the smallest cell in column major order 
    std::priority_queue<CellWrapper, std::vector<CellWrapper>, ColumnMajorCellCompareGT> m_cell_wrapper_pq;
#endif
};

#ifdef HTSDIR
class LoaderCombinedGVCFOperator : public LoaderOperatorBase
{
  public:
    LoaderCombinedGVCFOperator(const VidMapper* id_mapper, const std::string& config_filename,
        int partition_idx,
        const bool vid_mapper_file_required);
    virtual ~LoaderCombinedGVCFOperator()
    {
      clear();
      if(m_schema)
        delete m_schema;
      if(m_query_processor)
        delete m_query_processor;
      if(m_operator)
        delete m_operator;
      if(m_cell)
        delete m_cell;
      delete m_vcf_adapter;
    }
    virtual void operate(const void* cell_ptr);
    virtual bool overflow() const { return m_vcf_adapter->overflow(); }
    virtual void flush_output()
    {
      if(m_offload_vcf_output_processing)
        m_buffered_vcf_adapter->do_output();
    }
    virtual void post_operate_sequential()
    {
      if(m_offload_vcf_output_processing)
        m_buffered_vcf_adapter->advance_write_idx();
    }
    virtual void finish(const int64_t column_interval_end);
    void clear();
  private:
    VariantArraySchema* m_schema;
    VariantQueryProcessor* m_query_processor;
    const VidMapper* m_vid_mapper;
    VariantQueryConfig m_query_config;
    //Configuration for VCF adapter
    bool m_offload_vcf_output_processing;
    VCFAdapter* m_vcf_adapter;
    BufferedVCFAdapter* m_buffered_vcf_adapter; //points to same object as above
    //Operator that produces combined VCF record
    SingleVariantOperatorBase* m_operator;
    Variant m_variant;
    BufferVariantCell* m_cell;
    //PQ and aux structures
    VariantCallEndPQ m_end_pq;
    std::vector<VariantCall*> m_tmp_pq_vector;
    //Position trackers
    int64_t m_current_start_position;
    int64_t m_next_start_position;
    //Deletions
    uint64_t m_num_calls_with_deletions;
    //Profiling stat
    GTProfileStats m_stats;
    GTProfileStats* m_stats_ptr;
#ifdef DO_MEMORY_PROFILING
    size_t m_next_memory_limit;
#endif
};
#endif

#endif
