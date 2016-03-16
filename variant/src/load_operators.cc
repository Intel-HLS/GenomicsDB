#include "load_operators.h"
#include "json_config.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw LoadOperatorException(#X);

LoaderArrayWriter::LoaderArrayWriter(const VidMapper* id_mapper, const std::string& config_filename, int rank)
  : LoaderOperatorBase(), m_array_descriptor(-1), m_schema(0), m_storage_manager(0)
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
  bool row_based_partitioning = (json_doc.HasMember("row_based_partitioning") && json_doc["row_based_partitioning"].GetBool());
  RowRange row_partition(0, id_mapper->get_num_callsets()-1);
  if(row_based_partitioning)
  {
    row_partition = json_config.get_row_partition(rank);
    row_partition.second = std::min(row_partition.second, static_cast<int64_t>(id_mapper->get_num_callsets()-1));
  }
  id_mapper->build_tiledb_array_schema(m_schema, array_name, row_based_partitioning, row_partition, true);
  //Storage manager
  m_storage_manager = new VariantStorageManager(workspace);
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
#ifdef DUPLICATE_CELL_AT_END
  m_cell_copies.clear();
  m_last_end_position_for_row.resize(id_mapper->get_num_callsets(), -1ll);
#endif
}

#ifdef DUPLICATE_CELL_AT_END
void LoaderArrayWriter::write_top_element_to_disk()
{
  //Copy not reference
  CellWrapper top_element = m_cell_wrapper_pq.top();
  m_cell_wrapper_pq.pop();
  auto idx_in_vector = top_element.m_idx_in_cell_copies_vector;
  m_storage_manager->write_cell_sorted(m_array_descriptor,
      reinterpret_cast<const void*>(m_cell_copies[idx_in_vector]));
  //If this is a begin cell and spans multiple columns, retain this copy for the END in the PQ
  if(top_element.m_end_column > top_element.m_begin_column)
  {
    //swap begin/end
    std::swap<int64_t>(top_element.m_begin_column, top_element.m_end_column);
    //Update co-ordinate and END in the cell buffer
    auto* copy_ptr = m_cell_copies[idx_in_vector];
    //column is second co-ordinate
    *(reinterpret_cast<int64_t*>(copy_ptr+sizeof(int64_t))) = top_element.m_begin_column;
    //END is after co-ordinates and cell_size
    *(reinterpret_cast<int64_t*>(copy_ptr+2*sizeof(int64_t)+sizeof(size_t))) = top_element.m_end_column;
    //Add to PQ again
    m_cell_wrapper_pq.push(top_element);
  }
  else  //no need to keep this cell anymore, free "heap"
    m_memory_manager.push(idx_in_vector);
}
#endif

void LoaderArrayWriter::operate(const void* cell_ptr)
{
  assert(m_storage_manager);
#ifdef DUPLICATE_CELL_AT_END
  const uint8_t* ptr = reinterpret_cast<const uint8_t*>(cell_ptr);
  //row
  auto row = *(reinterpret_cast<const int64_t*>(ptr));
  //column is second co-ordinate
  auto column_begin = *(reinterpret_cast<const int64_t*>(ptr+sizeof(int64_t)));
  //END is after co-ordinates and cell_size
  auto column_end = *(reinterpret_cast<const int64_t*>(ptr+2*sizeof(int64_t)+sizeof(size_t)));
  assert(column_end >= column_begin && static_cast<size_t>(row) < m_last_end_position_for_row.size());
  //cell size is after co-ordinates
  auto cell_size = *(reinterpret_cast<const size_t*>(ptr+2*sizeof(int64_t)));
  //Reason: the whole setup works only if the intervals for a given row/sample are non-overlapping. This
  //property must be enforced by the loader
  //We maintain the last END value seen for every row - if the new cell has a begin
  // <= last END, then we truncate the END to column_begin-1 and write out the END copy cell
  //However, we must ensure that the copy cell with the truncated END respects column order
  //Hence, we only write to disks those cells (and END cell copies) which are less than (row+1, column_begin-1)
  //That way a truncated cell copy will be inserted at the correct position
  //Note that this increases memory consumption and run-time as every cell needs to be copied here
  CellWrapper curr_cell_wrapper({row+1, column_begin-1, -1, 0ull});
  ColumnMajorCellCompareGT cmp_op_GT;
  //Loop till (row+1, column_begin-1) > PQ top and write the top element to disk
  while(!m_cell_wrapper_pq.empty() && cmp_op_GT(curr_cell_wrapper, m_cell_wrapper_pq.top()))
    write_top_element_to_disk();
  //Check whether the last END value for this row overlaps current cell
  //If yes, must flush the entry from the PQ, update its co-ordinate to be column_begin-1 and insert into PQ again
  //Hopefully, entering this if statement is NOT the common case
  if(m_last_end_position_for_row[row] >= column_begin)
  {
    std::vector<CellWrapper> tmp_wrapper_vector;
    auto found_element = false;
    while(!m_cell_wrapper_pq.empty() && !found_element)
    {
      auto& top_ref = m_cell_wrapper_pq.top();
      if(top_ref.m_row == row)
        found_element = true;
      tmp_wrapper_vector.push_back(top_ref);
      m_cell_wrapper_pq.pop();
    }
    assert(found_element && tmp_wrapper_vector.size() > 0u);
    //Handle the last element separately
    for(auto i=0ull;i+1u<tmp_wrapper_vector.size();++i)
      m_cell_wrapper_pq.push(tmp_wrapper_vector[i]);
    //The cell corresponding to this row
    auto& last_element = tmp_wrapper_vector.back();
    auto idx_in_vector = last_element.m_idx_in_cell_copies_vector;
    auto copy_ptr = m_cell_copies[idx_in_vector];
    //Should always be an END copy cell - why? Because if this is a valid begin cell, then
    //m_begin_column < column_begin and the cell would have been written to disk by the loop over
    //the PQ above
    if(last_element.m_end_column < last_element.m_begin_column) //end copy, update co-ordinate
    {
      assert(last_element.m_begin_column == m_last_end_position_for_row[row]);
      last_element.m_begin_column = column_begin-1;
      //if the END copy still is at a column > its begin position, then write to disk
      //Due to the way the loop over the PQ operates above, this END copy cell is the next cell to go to disk
      //If the END copy cell is at column == its begin position, then after truncation, the copy has become a 
      //single position cell and there is no need to write the END copy
      if(last_element.m_begin_column != last_element.m_end_column)
      {
        //column is second co-ordinate
        *(reinterpret_cast<int64_t*>(copy_ptr+sizeof(int64_t))) = last_element.m_begin_column;
        m_storage_manager->write_cell_sorted(m_array_descriptor, reinterpret_cast<const void*>(copy_ptr));
      }
      m_memory_manager.push(idx_in_vector); //"free" memory
    }
    else      //m_begin_column>=m_end_column>=column_begin, incorrect input data
      throw LoadOperatorException(std::string("ERROR: two cells in incorrect order found\nPrevious cell: ")+
          std::to_string(last_element.m_row)+", "+std::to_string(last_element.m_begin_column)+", "+
          std::to_string(last_element.m_end_column)+
          "\nNew cell: "+std::to_string(row)+", "+std::to_string(column_begin));
  }
  size_t idx_in_vector = 0ull;
  if(m_memory_manager.empty())        //no free entries, need to allocate a new block
  {
    idx_in_vector = m_cell_copies.size();
    m_cell_copies.push_back(static_cast<uint8_t*>(malloc(cell_size)));
  }
  else        //free entry
  {
    idx_in_vector = m_memory_manager.top();
    m_memory_manager.pop();
    //realloc to new size
    m_cell_copies[idx_in_vector] = static_cast<uint8_t*>(realloc(m_cell_copies[idx_in_vector], cell_size));
  }
  VERIFY_OR_THROW(m_cell_copies[idx_in_vector] && "Memory allocation failed while creating copy of cell");
  memcpy(m_cell_copies[idx_in_vector], ptr, cell_size);
  //Update the cell wrapper structure
  curr_cell_wrapper.m_row = row;
  curr_cell_wrapper.m_begin_column = column_begin;
  curr_cell_wrapper.m_end_column = column_end;
  curr_cell_wrapper.m_idx_in_cell_copies_vector = idx_in_vector;
  //insert CellWrapper pointer into PQ
  m_cell_wrapper_pq.push(curr_cell_wrapper);
  //Update last END value seen
  m_last_end_position_for_row[row] = column_end;
#else //ifdef DUPLICATE_CELL_AT_END
  m_storage_manager->write_cell_sorted(m_array_descriptor, cell_ptr);
#endif //ifdef DUPLICATE_CELL_AT_END
}

void LoaderArrayWriter::finish(const int64_t column_interval_end)
{
#ifdef DUPLICATE_CELL_AT_END
  //some cells may be left in the PQ, write them to disk
  while(!m_cell_wrapper_pq.empty())
    write_top_element_to_disk();
#endif
  if(m_storage_manager && m_array_descriptor >= 0)
    m_storage_manager->close_array(m_array_descriptor);
}

#ifdef HTSDIR
LoaderCombinedGVCFOperator::LoaderCombinedGVCFOperator(const VidMapper* id_mapper, const std::string& config_filename,
    bool handle_spanning_deletions, int partition_idx, const ColumnRange& partition_range)
  : LoaderOperatorBase(), m_schema(0), m_query_processor(0), m_operator(0)
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
  m_vid_mapper->build_tiledb_array_schema(m_schema, "", false, RowRange(0, id_mapper->get_num_callsets()-1), false);
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
  m_cell = new BufferVariantCell(*m_schema, m_query_config);
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
  m_cell->set_cell(cell_ptr);
  m_query_processor->scan_handle_cell(m_query_config, 0u,
      m_variant, *m_operator, *m_cell,
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
