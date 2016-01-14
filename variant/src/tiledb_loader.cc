#include "tiledb_loader.h"
#include "rapidjson/document.h"
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw VCF2TileDBException(#X);

void LoaderConverterMessageExchange::resize_vectors(int num_divisions, int64_t total_size)
{
  m_all_num_tiledb_row_idx_vec_request.resize(num_divisions);
  m_all_num_tiledb_row_idx_vec_response.resize(num_divisions);
  for(auto i=0;i<num_divisions;++i)
    m_all_num_tiledb_row_idx_vec_request[i] = m_all_num_tiledb_row_idx_vec_response[i] = 0;
  m_all_tiledb_row_idx_vec_request.resize(total_size);
  m_all_tiledb_row_idx_vec_response.resize(total_size);
}

void LoaderConverterMessageExchange::initialize_from_converter(int num_partitions, int64_t num_owned_callsets)
{
  resize_vectors(num_partitions, num_partitions*num_owned_callsets);
  m_max_num_values_per_division.resize(num_partitions);
  m_idx_offset_per_division.resize(num_partitions);
  auto idx_offset = 0ull;
  for(auto i=0ull;i<m_max_num_values_per_division.size();++i)
  {
    m_max_num_values_per_division[i] = num_owned_callsets;
    m_idx_offset_per_division[i] = idx_offset;
    idx_offset += num_owned_callsets;
  }
}
   
//Same as converter - single partition
void LoaderConverterMessageExchange::initialize_from_loader(int64_t all_callsets)
{
  initialize_from_converter(1, all_callsets);
}

VCF2TileDBLoaderConverterBase::VCF2TileDBLoaderConverterBase(const std::string& config_filename, int idx)
{
  clear();
  m_idx = idx;
  //Parse json configuration
  rapidjson::Document json_doc;
  std::ifstream ifs(config_filename.c_str());
  VERIFY_OR_THROW(ifs.is_open());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  json_doc.Parse(str.c_str());
  //Column partitions
  VERIFY_OR_THROW(json_doc.HasMember("column_partitions"));
  {
    auto& column_partitions_dict = json_doc["column_partitions"];
    //column_partitions_dict itself a dictionary of the form { "0" : { "begin" : <value> } }
    VERIFY_OR_THROW(column_partitions_dict.IsObject());
    m_column_partition_begin_values.resize(column_partitions_dict.MemberCount());
    auto partition_idx = 0u;
    for(auto b=column_partitions_dict.MemberBegin(), e=column_partitions_dict.MemberEnd();b!=e;++b,++partition_idx)
    {
      const auto& curr_obj = *b;
      //{ "begin" : <Val> }
      const auto& curr_partition_info_dict = curr_obj.value;
      VERIFY_OR_THROW(curr_partition_info_dict.IsObject());
      VERIFY_OR_THROW(curr_partition_info_dict.HasMember("begin"));
      m_column_partition_begin_values[partition_idx] = curr_partition_info_dict["begin"].GetInt64();
    }
    //Sort in ascending order
    std::sort(m_column_partition_begin_values.begin(), m_column_partition_begin_values.end());
  }
  //Must have path to vid_mapping_file
  VERIFY_OR_THROW(json_doc.HasMember("vid_mapping_file"));
  m_vid_mapping_filename = json_doc["vid_mapping_file"].GetString();
  //Buffer size per column partition
  VERIFY_OR_THROW(json_doc.HasMember("size_per_column_partition"));
  m_per_partition_size = json_doc["size_per_column_partition"].GetInt64();
  //Obtain number of converters
  m_num_converter_processes = 0;
  if(json_doc.HasMember("num_converter_processes"))
    m_num_converter_processes = json_doc["num_converter_processes"].GetInt64();
  //Converter processes run independent of loader when num_converter_processes > 0
  m_standalone_converter_process = (m_num_converter_processes) ? true : false;
  //treat deletions as intervals
  if(json_doc.HasMember("treat_deletions_as_intervals"))
    m_treat_deletions_as_intervals = json_doc["treat_deletions_as_intervals"].GetBool();
  else
    m_treat_deletions_as_intervals = false;
  //Circular buffers - 3 needed in non-standalone converter mode
  resize_circular_buffers(3u);
  //Exchange structure
  m_owned_exchanges.resize(2u);
  //callset mapping file - check if defined in upper level json
  m_callset_mapping_file = "";
  if(json_doc.HasMember("callset_mapping_file") && json_doc["callset_mapping_file"].IsString())
    m_callset_mapping_file = json_doc["callset_mapping_file"].GetString();
  //Produce combined vcf
  m_produce_combined_vcf = false;
  if(json_doc.HasMember("produce_combined_vcf") && json_doc["produce_combined_vcf"].GetBool())
    m_produce_combined_vcf = true;
}

void VCF2TileDBLoaderConverterBase::clear()
{
  m_vid_mapping_filename.clear();
  m_callset_mapping_file.clear();
  m_column_partition_begin_values.clear();
  m_ping_pong_buffers.clear();
  m_owned_exchanges.clear();
}

VCF2TileDBConverter::VCF2TileDBConverter(const std::string& config_filename, int idx, VidMapper* vid_mapper, 
    std::vector<std::vector<uint8_t>>* buffers, std::vector<LoaderConverterMessageExchange>* exchange_vector)
  : VCF2TileDBLoaderConverterBase(config_filename, idx)
{
  m_vid_mapper = 0;
  clear();
  //Converter processes run independent of loader when num_converter_processes > 0
  if(m_standalone_converter_process)
  {
    VERIFY_OR_THROW(m_idx < m_num_converter_processes);
    //For standalone processes, must initialize VidMapper
    m_vid_mapper = static_cast<VidMapper*>(new FileBasedVidMapper(m_vid_mapping_filename, m_callset_mapping_file));
    //2 entries sufficient
    resize_circular_buffers(2u);
    //Buffer "pointers"
    m_cell_data_buffers.resize(m_num_entries_in_circular_buffer);
    for(auto i=0u;i<m_ping_pong_buffers.size();++i)
      m_cell_data_buffers[i] = &(m_ping_pong_buffers[i]);
    //Exchange pointers
    m_exchanges.resize(m_owned_exchanges.size());
    for(auto i=0u;i<m_owned_exchanges.size();++i)
      m_exchanges[i] = &(m_owned_exchanges[i]);
  }
  else
  {
    VERIFY_OR_THROW(static_cast<size_t>(m_idx) < m_column_partition_begin_values.size());
    m_vid_mapper = vid_mapper;
    VERIFY_OR_THROW(m_vid_mapper);
    //Buffer maintained external to converter, only maintain pointers
    VERIFY_OR_THROW(buffers);
    m_num_entries_in_circular_buffer = buffers->size();
    m_cell_data_buffers.resize(buffers->size());
    for(auto i=0u;i<m_cell_data_buffers.size();++i)
      m_cell_data_buffers[i] = &((*buffers)[i]);
    //Exchanges maintained external to converter, only maintain pointers
    VERIFY_OR_THROW(exchange_vector);
    m_exchanges.resize(exchange_vector->size());
    for(auto i=0u;i<m_exchanges.size();++i)
      m_exchanges[i] = &((*exchange_vector)[i]);
  } 
  m_max_size_per_callset = m_per_partition_size/m_vid_mapper->get_num_callsets();
  initialize_vcf2binary_objects();
  initialize_column_batch_objects();
  //For standalone converter objects, allocate ping-pong buffers and exchange objects
  if(m_standalone_converter_process)
  {
    for(auto& x : m_ping_pong_buffers)
      x.resize(m_max_size_per_callset*m_num_callsets_owned*m_partition_batch.size());
    for(auto& x : m_owned_exchanges)
      x.initialize_from_converter(m_partition_batch.size(), m_num_callsets_owned);
  }
}

VCF2TileDBConverter::~VCF2TileDBConverter()
{
  clear();
  if(m_standalone_converter_process && m_vid_mapper)
    delete m_vid_mapper;
  m_vid_mapper = 0;
}

void VCF2TileDBConverter::clear()
{
  m_cell_data_buffers.clear();
  m_partition_batch.clear();
  m_vcf_fields.clear();
  m_vcf2binary_handlers.clear();
  m_exchanges.clear();
}

void VCF2TileDBConverter::initialize_vcf2binary_objects()
{
  //VCF fields
  m_vid_mapper->build_vcf_fields_vectors(m_vcf_fields);
  //If standalone, deal only with subset of files assigned to this converter
  if(m_standalone_converter_process)
  {
    //Get list of files handled by this converter
    auto& global_file_idx_vec = m_vid_mapper->get_global_file_idxs_owned_by(m_idx);
    for(auto i=0ull;i<global_file_idx_vec.size();++i)
    {
      auto global_file_idx = global_file_idx_vec[i];
      auto& file_info = m_vid_mapper->get_file_info(global_file_idx);
      m_vcf2binary_handlers.emplace_back( 
          file_info.m_name, m_vcf_fields, i, *m_vid_mapper, m_column_partition_begin_values,
          m_max_size_per_callset,
          m_treat_deletions_as_intervals, false, false, false
          );
    }
  }
  else
  {
    //Same process as loader - must read all files
    //Also, only 1 partition needs to be handled  - the column partition corresponding to the loader
    auto partition_bounds=std::vector<int64_t>(2u);
    partition_bounds[0] = m_column_partition_begin_values[m_idx];
    if(static_cast<size_t>(m_idx) < m_column_partition_begin_values.size()-1)
      partition_bounds[1] = m_column_partition_begin_values[m_idx+1];
    else
      partition_bounds[1] = INT64_MAX;
    for(auto i=0ll;i<m_vid_mapper->get_num_files();++i)
    {
      auto global_file_idx = i;
      auto& file_info = m_vid_mapper->get_file_info(global_file_idx);
      m_vcf2binary_handlers.emplace_back(
          file_info.m_name, m_vcf_fields, i, *m_vid_mapper, partition_bounds,
          m_max_size_per_callset,
          m_treat_deletions_as_intervals, false, false, false
          );
    }
  }
}
  
void VCF2TileDBConverter::initialize_column_batch_objects()
{
  std::vector<int64_t> num_callsets_in_file;
  //If standalone, deal only with subset of files assigned to this converter
  if(m_standalone_converter_process)
  {
    //Get list of files handled by this converter
    auto& global_file_idx_vec = m_vid_mapper->get_global_file_idxs_owned_by(m_idx);
    num_callsets_in_file.resize(global_file_idx_vec.size());
    for(auto i=0ull;i<global_file_idx_vec.size();++i)
    {
      auto global_file_idx = global_file_idx_vec[i];
      auto& file_info = m_vid_mapper->get_file_info(global_file_idx);
      num_callsets_in_file[i] = file_info.get_num_callsets();
    }
  }
  else
  {
    //Same process as loader - must read all files
    num_callsets_in_file.resize(m_vid_mapper->get_num_files());
    for(auto i=0ll;i<m_vid_mapper->get_num_files();++i)
    {
      auto global_file_idx = i;
      auto& file_info = m_vid_mapper->get_file_info(global_file_idx);
      num_callsets_in_file[i] = file_info.get_num_callsets();
    }
  }
  //If standalone converter process, initialize for all partitions
  //Else only allocate single column partition idx corresponding to the loader
  auto num_column_partitions = m_standalone_converter_process ? m_column_partition_begin_values.size() : 1u;
  for(auto i=0u;i<num_column_partitions;++i)
    m_partition_batch.emplace_back(i, m_max_size_per_callset, num_callsets_in_file, m_num_entries_in_circular_buffer);
  m_num_callsets_owned = 0;
  for(auto x : num_callsets_in_file)
    m_num_callsets_owned += x;
  //Update row idx to ordering
  m_tiledb_row_idx_to_order = std::move(std::vector<int64_t>(m_num_callsets_owned, -1ll));
  int64_t order = 0ll;
  for(const auto& x : m_vcf2binary_handlers)
    x.set_order_of_enabled_callsets(order, m_tiledb_row_idx_to_order);
  assert(order == m_num_callsets_owned);
}

void VCF2TileDBConverter::activate_next_batch(const unsigned exchange_idx, const int partition_idx)
{
  assert(exchange_idx < m_exchanges.size());
  auto& curr_exchange = *(m_exchanges[exchange_idx]);
  auto& all_partitions_tiledb_row_idx_vec = curr_exchange.m_all_tiledb_row_idx_vec_request;
  assert(static_cast<size_t>(partition_idx) < m_partition_batch.size());
  auto idx_offset = curr_exchange.get_idx_offset_for_partition(partition_idx);
  for(auto i=0ll;i<curr_exchange.m_all_num_tiledb_row_idx_vec_request[partition_idx];++i)
  {
    assert(static_cast<size_t>(idx_offset+i) < all_partitions_tiledb_row_idx_vec.size());
    auto row_idx = all_partitions_tiledb_row_idx_vec[idx_offset+i];
    int64_t local_file_idx = -1;
    m_vid_mapper->get_local_file_idx_for_row(row_idx, local_file_idx);
    assert(local_file_idx >= 0 && static_cast<size_t>(local_file_idx) < m_vcf2binary_handlers.size());
    //Advance write idx
    //If non-standalone converter, advance read idx also
    m_partition_batch[partition_idx].activate_file(local_file_idx, !m_standalone_converter_process);
  }
  //Compute offsets in buffer if standalone, otherwise maintain offsets computed during initialization
  if(m_standalone_converter_process)
    m_partition_batch[partition_idx].update_buffer_offsets();
}

void VCF2TileDBConverter::read_next_batch(const unsigned exchange_idx)
{
  assert(exchange_idx < m_exchanges.size());
  auto& curr_exchange = *(m_exchanges[exchange_idx]);
  for(auto partition_idx=0u;partition_idx<m_partition_batch.size();++partition_idx)
    if(curr_exchange.is_partition_requested_by_loader(partition_idx))
    {
      activate_next_batch(exchange_idx, partition_idx);
      //Find callsets which still have new data in this partition
      int64_t idx_offset = curr_exchange.get_idx_offset_for_partition(partition_idx);
      auto begin_idx_offset = idx_offset;
      for(auto& vcf_handler : m_vcf2binary_handlers)
        vcf_handler.list_active_row_idxs(m_partition_batch[partition_idx], idx_offset,
            curr_exchange.m_all_tiledb_row_idx_vec_response);
      curr_exchange.m_all_num_tiledb_row_idx_vec_response[partition_idx] = idx_offset - begin_idx_offset;
    }
#pragma omp parallel for default(shared)
  for(auto i=0u;i<m_vcf2binary_handlers.size();++i)
    m_vcf2binary_handlers[i].read_next_batch(m_cell_data_buffers, m_partition_batch, false);
  curr_exchange.m_is_serviced = true;
}

void VCF2TileDBConverter::dump_latest_buffer(unsigned exchange_idx, std::ostream& osptr) const
{
  auto& curr_exchange = *(m_exchanges[exchange_idx]);
  osptr << "Batch in exchange "<<exchange_idx<<"\n";
  for(auto partition_idx=0u;partition_idx<m_partition_batch.size();++partition_idx)
  {
    if(curr_exchange.is_partition_requested_by_loader(partition_idx))
    {
      int64_t idx_offset = curr_exchange.get_idx_offset_for_partition(partition_idx);
      for(auto i=0ll;i<curr_exchange.m_all_num_tiledb_row_idx_vec_response[partition_idx];++i)
      {
        auto row_idx = curr_exchange.m_all_tiledb_row_idx_vec_response[idx_offset+i];
        int64_t local_file_idx = -1;
        auto status = 
          m_vid_mapper->get_local_file_idx_for_row(row_idx, local_file_idx);
        assert(status);
        auto& file_batch = m_partition_batch[partition_idx].get_partition_file_batch(local_file_idx);
        auto order = m_standalone_converter_process ? i : m_tiledb_row_idx_to_order[row_idx];
        assert(order >= 0);
        auto offset = m_partition_batch[partition_idx].get_partition_begin_offset() + order*m_max_size_per_callset;
        for(auto j=0ll;j<m_max_size_per_callset;++j)
        {
          auto val = (*(m_cell_data_buffers[file_batch.get_buffer_idx()]))[offset+j];
          if(val != 0)
            osptr << static_cast<char>(val);
          else
            break;
        }
      }
    }
  }
}

//Loader functions
VCF2TileDBLoader::VCF2TileDBLoader(const std::string& config_filename, int idx)
  : VCF2TileDBLoaderConverterBase(config_filename, idx)
{
  m_converter = 0;
  clear();
  VERIFY_OR_THROW(static_cast<size_t>(m_idx) < m_column_partition_begin_values.size());
  m_vid_mapper = static_cast<VidMapper*>(new FileBasedVidMapper(m_vid_mapping_filename, m_callset_mapping_file));
  m_max_size_per_callset = m_per_partition_size/m_vid_mapper->get_num_callsets();
  //Converter processes run independent of loader when num_converter_processes > 0
  if(m_standalone_converter_process)
  {
    resize_circular_buffers(4u);
    //Allocate exchange objects
  }
  else
  {
    //Allocate exchange objects
    for(auto& x : m_owned_exchanges)
    {
      x.initialize_from_loader(m_vid_mapper->get_num_callsets());
      x.m_all_num_tiledb_row_idx_vec_request[0] = m_vid_mapper->get_num_callsets();
      for(auto i=0ll;i<m_vid_mapper->get_num_callsets();++i)
        x.m_all_tiledb_row_idx_vec_request[x.get_idx_offset_for_converter(0)+i] = i;
    }
    m_converter = new VCF2TileDBConverter(config_filename, idx, m_vid_mapper, &m_ping_pong_buffers, &m_owned_exchanges);
  }
  //Allocate buffers
  for(auto i=0u;i<m_ping_pong_buffers.size();++i)
    m_ping_pong_buffers[i].resize(m_per_partition_size);
  //Circular buffer control
  m_row_idx_to_buffer_control.resize(m_vid_mapper->get_num_callsets(), CircularBufferController(m_num_entries_in_circular_buffer));
  //Priority queue elements
  m_pq_vector.resize(m_vid_mapper->get_num_callsets());
  m_rows_not_in_pq.resize(m_vid_mapper->get_num_callsets());
  for(auto i=0ull;i<m_pq_vector.size();++i)
  {
    m_pq_vector[i].m_row_idx = i;
    m_pq_vector[i].m_offset = get_buffer_start_offset_for_row_idx(i);
    m_rows_not_in_pq[i] = i;
  }
  //Operators
  if(m_produce_combined_vcf)
    m_operators.push_back(dynamic_cast<LoaderOperatorBase*>(
          new LoaderCombinedGVCFOperator(m_vid_mapper, config_filename, m_treat_deletions_as_intervals)));
}

#include "omp.h"
void VCF2TileDBLoader::read_all()
{
  auto num_exchanges = m_owned_exchanges.size();
  auto exchange_counter = num_exchanges-1u;
  while(true)
  {
    auto done=false;
    auto fetch_exchange_counter = (exchange_counter+1u)%num_exchanges;
    auto load_exchange_counter = exchange_counter;
    //For row idx requested, reserve entries
    reserve_entries_in_circular_buffer(fetch_exchange_counter);
#pragma omp parallel sections default(shared)
    {
#pragma omp section
      {
        m_converter->read_next_batch(fetch_exchange_counter);
      }
#pragma omp section
      {
        //done = dump_latest_buffer(load_exchange_counter, std::cout);
        done = produce_cells_in_column_major_order(load_exchange_counter);
      }
    }
    advance_write_idxs(fetch_exchange_counter);
    if(done)
      break;
    exchange_counter = (exchange_counter+1u)%num_exchanges;
  }
  for(auto op : m_operators)
    op->finish(get_column_partition_end());
}

void VCF2TileDBLoader::reserve_entries_in_circular_buffer(unsigned exchange_idx)
{
  auto converter_idx = 0u;
  auto& curr_exchange = m_owned_exchanges[exchange_idx];
  //Reserve entry in circular buffer - implies that entry is reserved, but no valid data exists
  int64_t idx_offset = curr_exchange.get_idx_offset_for_converter(converter_idx);
  for(auto i=0ll;i<curr_exchange.m_all_num_tiledb_row_idx_vec_request[converter_idx];++i)
  {
    auto row_idx = curr_exchange.m_all_tiledb_row_idx_vec_request[idx_offset+i];
    assert(m_row_idx_to_buffer_control[row_idx].get_num_empty_entries() > 0u);
    m_row_idx_to_buffer_control[row_idx].reserve_entry();
  }
}

void VCF2TileDBLoader::advance_write_idxs(unsigned exchange_idx)
{
  auto converter_idx = 0u;
  auto& curr_exchange = m_owned_exchanges[exchange_idx];
  if(!curr_exchange.m_is_serviced)
    return;
  //Advance circular buffer control
  int64_t idx_offset = curr_exchange.get_idx_offset_for_converter(converter_idx);
  for(auto i=0ll;i<curr_exchange.m_all_num_tiledb_row_idx_vec_response[converter_idx];++i)
  {
    auto row_idx = curr_exchange.m_all_tiledb_row_idx_vec_response[idx_offset+i];
    //Advance and un-reserve
    m_row_idx_to_buffer_control[row_idx].advance_write_idx(true);
  }
}

bool VCF2TileDBLoader::dump_latest_buffer(unsigned exchange_idx, std::ostream& osptr)
{
  auto& curr_exchange = m_owned_exchanges[exchange_idx];
  auto converter_idx = 0u;
  if(curr_exchange.m_is_serviced)
  {
    osptr << "Batch in exchange "<<exchange_idx<<"\n";
    int64_t idx_offset = curr_exchange.get_idx_offset_for_converter(converter_idx);
    for(auto i=0ll;i<curr_exchange.m_all_num_tiledb_row_idx_vec_response[converter_idx];++i)
    {
      auto row_idx = curr_exchange.m_all_tiledb_row_idx_vec_response[idx_offset+i];
      assert(row_idx >= 0 && static_cast<size_t>(row_idx) < m_row_idx_to_buffer_control.size());
      //Add to next request
      curr_exchange.m_all_tiledb_row_idx_vec_request[idx_offset + i] = row_idx;
      //Ping pong buffering control
      auto buffer_idx = m_row_idx_to_buffer_control[row_idx].get_read_idx();
      assert(m_row_idx_to_buffer_control[row_idx].get_num_entries_with_valid_data() > 0u);
      auto order = get_order_for_row_idx(row_idx);
      assert(order >= 0);
      int64_t offset = order*m_max_size_per_callset;
      auto j=0ll;
      for(;j<m_max_size_per_callset;++j)
      {
        auto val = m_ping_pong_buffers[buffer_idx][offset+j];
        if(val != 0)
          osptr << static_cast<char>(val);
        else
          break;
      }
      m_row_idx_to_buffer_control[row_idx].advance_read_idx();
    }
    curr_exchange.m_all_num_tiledb_row_idx_vec_request[converter_idx] = curr_exchange.m_all_num_tiledb_row_idx_vec_response[converter_idx];
    return (curr_exchange.m_all_num_tiledb_row_idx_vec_response[converter_idx] == 0);
  }
  else
    return false;
}

bool VCF2TileDBLoader::read_cell_from_buffer(const int64_t row_idx)
{
  assert(row_idx >= 0 && static_cast<size_t>(row_idx) < m_row_idx_to_buffer_control.size());
  auto& buffer_control = m_row_idx_to_buffer_control[row_idx];
  //No valid entries exist
  if(buffer_control.get_num_entries_with_valid_data() == 0u)
    return false;
  auto& pq_element = m_pq_vector[row_idx];
  //Past the buffer limit for this callset
  if(static_cast<size_t>(pq_element.m_offset) + sizeof(int64_t) > 
      static_cast<size_t>(get_buffer_start_offset_for_row_idx(row_idx)) + m_max_size_per_callset)
    return false;
  auto& curr_buffer = m_ping_pong_buffers[buffer_control.get_read_idx()];
  auto ptr = reinterpret_cast<const int64_t*>(&(curr_buffer[pq_element.m_offset]));
  auto row_idx_in_buffer = ptr[0];
  //row idx in buffer == NULL, no more valid data
  if(row_idx_in_buffer == get_tiledb_null_value<int64_t>())
    return false;
  assert(row_idx_in_buffer == row_idx);
  pq_element.m_column = ptr[1];
  return true;
}

bool VCF2TileDBLoader::read_next_cell_from_buffer(const int64_t row_idx)
{
  //curr position is valid
  assert(read_cell_from_buffer(row_idx));
  auto& buffer_control = m_row_idx_to_buffer_control[row_idx];
  auto& pq_element = m_pq_vector[row_idx];
  auto& curr_buffer = m_ping_pong_buffers[buffer_control.get_read_idx()];
  //Cell size is after the coordinates
  auto ptr = reinterpret_cast<const size_t*>(&(curr_buffer[pq_element.m_offset+2*sizeof(int64_t)]));
  auto cell_size = *ptr;
  //Update offset
  pq_element.m_offset += cell_size;
  //Check if next valid cell is within the buffer limit
  if(read_cell_from_buffer(row_idx))
    return true;
  //No valid cell found in current buffer - advance read idx
  buffer_control.advance_read_idx();
  //Reset offset
  pq_element.m_offset = get_buffer_start_offset_for_row_idx(row_idx);
  //If already crossed buffer once in this batch, don't bother reading further
  if(pq_element.m_crossed_one_buffer)
    return false;
  pq_element.m_crossed_one_buffer = true;
  return read_cell_from_buffer(row_idx);
}

bool VCF2TileDBLoader::produce_cells_in_column_major_order(unsigned exchange_idx)
{
  auto& curr_exchange = m_owned_exchanges[exchange_idx];
  if(!curr_exchange.m_is_serviced)
    return false;
  auto converter_idx = 0u;
  auto idx_offset = curr_exchange.get_idx_offset_for_converter(converter_idx);
  //Add callsets that are not in PQ into the PQ if valid cells found
  for(auto i=0ull;i<m_rows_not_in_pq.size();++i)
  {
    auto row_idx = m_rows_not_in_pq[i];
    auto valid_cell_found = read_cell_from_buffer(row_idx);
    if(valid_cell_found)
      m_column_major_pq.push(&(m_pq_vector[row_idx]));
    else
      m_pq_vector[row_idx].m_completed = true;
  }
  auto num_rows_not_in_pq = 0ull;
  //No re-allocation as resize() doesn't reduce capacity
  m_rows_not_in_pq.resize(m_vid_mapper->get_num_callsets());
  auto top_column = -1ll;
  auto hit_invalid_cell = false;
  while(!m_column_major_pq.empty() && (!hit_invalid_cell || (m_column_major_pq.top())->m_column == top_column))
  {
    auto* top_ptr = m_column_major_pq.top();
    m_column_major_pq.pop();
    auto row_idx = top_ptr->m_row_idx;
    //std::cout << row_idx <<","<<top_ptr->m_column<<"\n";
    const auto& buffer = m_ping_pong_buffers[m_row_idx_to_buffer_control[row_idx].get_read_idx()];
    auto offset = top_ptr->m_offset;
    for(auto op : m_operators)
      op->operate(reinterpret_cast<const void*>(&(buffer[offset])));
    auto valid_cell_found = read_next_cell_from_buffer(row_idx);
    if(valid_cell_found)
      m_column_major_pq.push(&(m_pq_vector[row_idx]));
    else
    {
      if(!hit_invalid_cell)     //first invalid cell found
      {
        hit_invalid_cell = true;
        top_column = top_ptr->m_column;
      }
      m_rows_not_in_pq[num_rows_not_in_pq++] = row_idx;
    }
  }
  //No re-allocation as resize() doesn't reduce capacity
  m_rows_not_in_pq.resize(num_rows_not_in_pq);
  //Find rows to request in next batch - rows which have empty space
  auto num_rows_in_next_request = 0ull;
  for(auto i=0ll;i<m_vid_mapper->get_num_callsets();++i)
  {
    auto& pq_element = m_pq_vector[i];
    pq_element.m_crossed_one_buffer = false;
    //Space in circular buffer
    if(!pq_element.m_completed && m_row_idx_to_buffer_control[i].get_num_empty_entries() > 0u)
    {
      curr_exchange.m_all_tiledb_row_idx_vec_request[idx_offset + num_rows_in_next_request] = i;
      ++num_rows_in_next_request;
    }
  }
  curr_exchange.m_all_num_tiledb_row_idx_vec_request[converter_idx] = num_rows_in_next_request;
  return (m_column_major_pq.empty() && num_rows_in_next_request == 0u && num_rows_not_in_pq == 0u);
}

void VCF2TileDBLoader::clear()
{
  m_row_idx_to_buffer_control.clear();
  m_pq_vector.clear();
  m_operators.clear();
}
