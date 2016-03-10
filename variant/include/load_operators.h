#ifndef LOAD_OPERATORS_H
#define LOAD_OPERATORS_H

#include "vid_mapper.h"
#include "query_variants.h"
#include "broad_combined_gvcf.h" 

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
    LoaderOperatorBase() { ; }
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
    virtual void finish(const int64_t column_interval_end) { ; }
};

class LoaderArrayWriter : public LoaderOperatorBase
{
  public:
    LoaderArrayWriter(const VidMapper* id_mapper, const std::string& config_filename, int rank);
    virtual ~LoaderArrayWriter()
    {
      if(m_schema)
        delete m_schema;
      if(m_storage_manager)
        delete m_storage_manager;
#ifdef DUPLICATE_CELL_AT_END
      for(auto ptr : m_end_point_cell_copies)
        free(ptr);
      m_end_point_cell_copies.clear();
#endif
    }
    virtual void operate(const void* cell_ptr);
    virtual void finish(const int64_t column_interval_end);
  private:
    int m_array_descriptor;
    ArraySchema* m_schema;
    StorageManager* m_storage_manager;
#ifdef DUPLICATE_CELL_AT_END
    //Copy of cells to store at the end position
    std::vector<uint8_t*> m_end_point_cell_copies;
    //Mimics behavior of program heap for cell copies - minimizes number of frees/reallocations
    std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> m_memory_manager;
    //For use in priority queue
    typedef struct
    {
      int64_t m_row;
      int64_t m_column;
      size_t m_idx_in_cell_copies_vector;
    } CellWrapper;
    struct ColumnMajorCellCompareGT 
    {
      bool operator()(const CellWrapper& a, const CellWrapper& b) const
      {
        return ((a.m_column > b.m_column) || (a.m_column == b.m_column && a.m_row > b.m_row));
      }
    };
    //top() points to Cell in m_end_point_cell_wrapper_vector which contains the cell with the smallest end point
    std::priority_queue<CellWrapper, std::vector<CellWrapper>, ColumnMajorCellCompareGT> m_end_point_cell_wrapper_pq;
#endif
};

#ifdef HTSDIR
class LoaderCombinedGVCFOperator : public LoaderOperatorBase
{
  public:
    LoaderCombinedGVCFOperator(const VidMapper* id_mapper, const std::string& config_filename, bool handle_spanning_deletions,
        int partition_idx, const ColumnRange& partition_range);
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
    ArraySchema* m_schema;
    VariantQueryProcessor* m_query_processor;
    const VidMapper* m_vid_mapper;
    VariantQueryConfig m_query_config;
    //Configuration for VCF adapter
    bool m_offload_vcf_output_processing;
    VCFAdapter* m_vcf_adapter;
    BufferedVCFAdapter* m_buffered_vcf_adapter; //points to same object as above
    //Operator that produces combined VCF record
    BroadCombinedGVCFOperator* m_operator;
    Variant m_variant;
    BufferVariantCell* m_cell;
    //Column interval bounds
    ColumnRange m_partition;
    //PQ and aux structures
    VariantCallEndPQ m_end_pq;
    std::vector<VariantCall*> m_tmp_pq_vector;
    //Position trackers
    int64_t m_current_start_position;
    int64_t m_next_start_position;
    //Deletions
    bool m_handle_spanning_deletions;
    uint64_t m_num_calls_with_deletions;
    //Profiling stat
    GTProfileStats m_stats;
    GTProfileStats* m_stats_ptr;
};
#endif

#endif
