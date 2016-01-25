#include "load_operators.h"
#include "json_config.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw LoadOperatorException(#X);

LoaderArrayWriter::LoaderArrayWriter(const VidMapper* id_mapper, const std::string& config_filename, int rank)
  : LoaderOperatorBase(), m_schema(0), m_storage_manager(0), m_array_descriptor(-1)
{
  //Parse json configuration
  rapidjson::Document json_doc;
  std::ifstream ifs(config_filename.c_str());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  json_doc.Parse(str.c_str());
  //recreate array flag
  bool recreate_array = (json_doc.HasMember("delete_and_create_tiledb_array")
      && json_doc["delete_and_create_tiledb_array"].GetBool()) ? true : false;
  JSONConfigBase json_config;
  json_config.read_from_file(config_filename);
  auto workspace = json_config.get_workspace(rank);
  auto array_name = json_config.get_array_name(rank);
  //Schema
  id_mapper->build_tiledb_array_schema(m_schema, true, array_name);
  //Storage manager
  m_storage_manager = new StorageManager(workspace);
  auto mode = recreate_array ? "w" : "a";
  //Check if array already exists
  m_array_descriptor = m_storage_manager->open_array(array_name, mode);
  //Array does not exist - define it first
  if(m_array_descriptor < 0)
  {
    VERIFY_OR_THROW(m_storage_manager->define_array(m_schema) == TILEDB_OK
        && "Could not define TileDB array");
    //Open array in write mode
    m_array_descriptor = m_storage_manager->open_array(array_name, "w");
  }
  VERIFY_OR_THROW(m_array_descriptor != -1 && "Could not open TileDB array for loading");
}

#ifdef HTSDIR
LoaderCombinedGVCFOperator::LoaderCombinedGVCFOperator(const VidMapper* id_mapper, const std::string& config_filename,
    bool handle_spanning_deletions, int partition_idx, const ColumnRange& partition_range)
  : LoaderOperatorBase(), m_operator(0), m_query_processor(0), m_schema(0)
{
  clear();
  //Parse json configuration
  rapidjson::Document json_doc;
  std::ifstream ifs(config_filename.c_str());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  json_doc.Parse(str.c_str());
  //initialize arguments
  m_vid_mapper = id_mapper;
  //initialize query processor
  m_vid_mapper->build_tiledb_array_schema(m_schema, false);
  m_query_processor = new VariantQueryProcessor(*m_schema);
  //Initialize query config
  std::vector<std::string> query_attributes(m_schema->attribute_num());
  for(auto i=0;i<m_schema->attribute_num();++i)
    query_attributes[i] = m_schema->attribute_name(i);
  m_query_config.set_attributes_to_query(query_attributes);
  m_query_processor->do_query_bookkeeping(*m_schema, m_query_config);
  //Initialize VCF adapter
  if(json_doc.HasMember("offload_vcf_output_processing") && json_doc["offload_vcf_output_processing"].GetBool())
  {
    m_offload_vcf_output_processing = true;
    //2 entries in circular buffer, max #entries to use in each line_buffer
    m_buffered_vcf_adapter = new BufferedVCFAdapter(2u, m_vid_mapper->get_num_callsets());
    m_vcf_adapter = dynamic_cast<VCFAdapter*>(m_buffered_vcf_adapter);
  }
  else
  {
    m_offload_vcf_output_processing = false;
    m_vcf_adapter = new VCFAdapter();
    m_buffered_vcf_adapter = 0;
  }
  JSONVCFAdapterConfig vcf_adapter_config;
  vcf_adapter_config.read_from_file(config_filename, *m_vcf_adapter, "", partition_idx);
  //Initialize operator
  m_operator = new BroadCombinedGVCFOperator(*m_vcf_adapter, *m_vid_mapper, m_query_config);
  //Initialize variant
  m_variant = std::move(Variant(&m_query_config));
  m_variant.resize_based_on_query();
  //Cell
  m_cell = new Cell(m_schema, m_query_config.get_query_attributes_schema_idxs(), 0, true);
  //Partition bounds
  m_partition = partition_range;
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
  m_variant.clear();
  m_tmp_pq_vector.clear();
}

void LoaderCombinedGVCFOperator::operate(const void* cell_ptr)
{
  auto coords = reinterpret_cast<const int64_t*>(cell_ptr);
  auto column_value = coords[1];
#ifndef NDEBUG
  auto ptr = reinterpret_cast<const uint8_t*>(cell_ptr);
  //END value is after cooords and cell_size
  ptr += 2*sizeof(int64_t)+sizeof(size_t);
  auto end_value = *(reinterpret_cast<const int64_t*>(ptr));
  //Must cross partition bound
  assert(end_value >= m_partition.first && column_value <= m_partition.second);
#endif
  //Either un-initialized or VariantCall interval starts before/at partition begin value
  if(m_current_start_position < 0 || column_value <= m_partition.first)
  {
    m_current_start_position = column_value;
    m_variant.set_column_interval(column_value, column_value);
  }
  else  //column_value > m_partition.first, check if m_current_start_position < m_partition.first
    if(m_current_start_position < m_partition.first)
    {
      m_current_start_position = m_partition.first;
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
  assert(!m_offload_vcf_output_processing || m_buffered_vcf_adapter->get_num_entries_with_valid_data() == 0u);
  pre_operate_sequential();
  //Fix start and next_start positions if necessary
  if(m_current_start_position < m_partition.first)
  {
    m_current_start_position = m_partition.first;
    m_variant.set_column_interval(m_current_start_position, m_current_start_position);
  }
  m_next_start_position = (column_interval_end == INT64_MAX) ? INT64_MAX : column_interval_end+1;
  m_query_processor->handle_gvcf_ranges(m_end_pq, m_query_config, m_variant, *m_operator,
      m_current_start_position, m_next_start_position, column_interval_end == INT64_MAX, m_num_calls_with_deletions);
  post_operate_sequential();
  flush_output();
}

#endif
