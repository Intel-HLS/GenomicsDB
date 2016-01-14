#include "load_operators.h"
#include "json_config.h"

LoaderCombinedGVCFOperator::LoaderCombinedGVCFOperator(const VidMapper* id_mapper, const std::string& config_filename,
    bool handle_spanning_deletions)
  : LoaderOperatorBase(), m_operator(0), m_query_processor(0), m_schema(0)
{
  clear();
  m_vid_mapper = id_mapper;
  //initialize query processor
  m_vid_mapper->build_tiledb_array_schema(m_schema);
  m_query_processor = new VariantQueryProcessor(*m_schema);
  //Initialize query config
  std::vector<std::string> query_attributes(m_schema->attribute_num());
  for(auto i=0;i<m_schema->attribute_num();++i)
    query_attributes[i] = m_schema->attribute_name(i);
  m_query_config.set_attributes_to_query(query_attributes);
  m_query_processor->do_query_bookkeeping(*m_schema, m_query_config);
  //Initialize VCF adapter
  JSONVCFAdapterConfig vcf_adapter_config;
  vcf_adapter_config.read_from_file(config_filename, m_vcf_adapter);
  //Initialize operator
  m_operator = new BroadCombinedGVCFOperator(m_vcf_adapter, *m_vid_mapper, m_query_config);
  //Initialize variant
  m_variant = std::move(Variant(&m_query_config));
  m_variant.resize_based_on_query();
  //Cell
  m_cell = new Cell(m_schema, m_query_config.get_query_attributes_schema_idxs(), 0, true);
  //PQ elements
  m_tmp_pq_vector.resize(m_query_config.get_num_rows_to_query());
  //Position elements
  m_current_start_position = -1ll;
  m_next_start_position = -1ll;
  //Deletion flags
  m_num_calls_with_deletions = 0;
  m_handle_spanning_deletions = handle_spanning_deletions;
  //Profiling
#ifdef DO_PROFILING
  m_stats_ptr = &m_stats;
#else
  m_stats_ptr = 0;
#endif
}

void LoaderCombinedGVCFOperator::clear()
{
  m_query_config.clear();
  m_vcf_adapter.clear();
  m_variant.clear();
  m_tmp_pq_vector.clear();
}

void LoaderCombinedGVCFOperator::operate(const void* cell_ptr)
{
  if(m_current_start_position < 0)      //un-initialized
  {
    auto coords = reinterpret_cast<const int64_t*>(cell_ptr);
    m_current_start_position = coords[1];
    m_variant.set_column_interval(m_current_start_position, m_current_start_position);
  }
  m_query_processor->scan_handle_cell(m_query_config, 0u,
      m_variant, *m_operator, *m_cell, cell_ptr,
      m_end_pq, m_tmp_pq_vector,
      m_current_start_position, m_next_start_position,
      m_num_calls_with_deletions, m_handle_spanning_deletions,
      m_stats_ptr);
  return;
}

void LoaderCombinedGVCFOperator::finish(const int64_t column_interval_end)
{
  m_next_start_position = (column_interval_end == INT64_MAX) ? INT64_MAX : column_interval_end+1;
  m_query_processor->handle_gvcf_ranges(m_end_pq, m_query_config, m_variant, *m_operator,
      m_current_start_position, m_next_start_position, column_interval_end == INT64_MAX, m_num_calls_with_deletions);
}
