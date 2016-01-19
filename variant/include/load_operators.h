#ifndef LOAD_OPERATORS_H
#define LOAD_OPERATORS_H

#include "vid_mapper.h"
#include "query_variants.h"
#include "broad_combined_gvcf.h" 

class LoaderOperatorBase
{
  public:
    LoaderOperatorBase() { ; }
    virtual ~LoaderOperatorBase() { ; }
    /*
     * Virtual function that must be overridden by sub-classes
     */
    virtual void operate(const void* cell_ptr) = 0;
    /*
     * Called at the end - the argument is the column interval end limit
     */
    virtual void finish(const int64_t column_interval_end) { ; }
};

#ifdef HTSDIR
class LoaderCombinedGVCFOperator : public LoaderOperatorBase
{
  public:
    LoaderCombinedGVCFOperator(const VidMapper* id_mapper, const std::string& config_filename, bool handle_spanning_deletions);
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
    }
    virtual void operate(const void* cell_ptr);
    virtual void finish(const int64_t column_interval_end);
    void clear();
  private:
    ArraySchema* m_schema;
    VariantQueryProcessor* m_query_processor;
    const VidMapper* m_vid_mapper;
    VariantQueryConfig m_query_config;
    VCFAdapter m_vcf_adapter;
    BroadCombinedGVCFOperator* m_operator;
    Variant m_variant;
    Cell* m_cell;
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
