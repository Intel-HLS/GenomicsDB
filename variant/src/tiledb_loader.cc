#include "tiledb_loader.h"
#include "rapidjson/document.h"
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw VCF2TileDBException(#X);

#define NUM_PING_PONG_BUFFERS 2u

void LoaderConverterMessageExchange::resize_vectors(int num_divisions, int64_t total_size)
{
  m_all_num_tiledb_row_idx_vec_request.resize(num_divisions);
  m_all_num_tiledb_row_idx_vec_response.resize(num_divisions);
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
  //Ping-pong buffers
  m_ping_pong_buffers.resize(NUM_PING_PONG_BUFFERS);
  //Exchanges
  m_owned_exchanges.resize(NUM_PING_PONG_BUFFERS);
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
}

void VCF2TileDBLoaderConverterBase::clear()
{
  m_column_partition_begin_values.clear();
  m_ping_pong_buffers.clear();
  m_owned_exchanges.clear();
  m_vid_mapping_filename.clear();
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
    m_vid_mapper = static_cast<VidMapper*>(new FileBasedVidMapper(m_vid_mapping_filename));
    //Buffer "pointers"
    m_cell_data_buffers.resize(NUM_PING_PONG_BUFFERS);
    //Exchange pointers
    m_exchanges.resize(NUM_PING_PONG_BUFFERS);
    for(auto i=0u;i<m_ping_pong_buffers.size();++i)
    {
      m_cell_data_buffers[i] = &(m_ping_pong_buffers[i]);
      m_exchanges[i] = &(m_owned_exchanges[i]);
    } 
  }
  else
  {
    VERIFY_OR_THROW(static_cast<size_t>(m_idx) < m_column_partition_begin_values.size());
    m_vid_mapper = vid_mapper;
    VERIFY_OR_THROW(m_vid_mapper);
    //Buffer maintained external to converter, only maintain pointers
    VERIFY_OR_THROW(buffers);
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
    m_partition_batch.emplace_back(i, m_max_size_per_callset, num_callsets_in_file);
  m_num_callsets_owned = 0;
  for(auto x : num_callsets_in_file)
    m_num_callsets_owned += x;
}

void VCF2TileDBConverter::activate_next_batch(const unsigned exchange_idx, const int partition_idx)
{
  assert(exchange_idx < m_exchanges.size());
  auto& curr_exchange = *(m_exchanges[exchange_idx]);
  auto& all_partitions_tiledb_row_idx_vec = curr_exchange.m_all_tiledb_row_idx_vec_request;
  assert(static_cast<size_t>(partition_idx) < m_partition_batch.size());
  int64_t begin_idx = curr_exchange.get_idx_offset_for_partition(partition_idx);
  for(auto i=0ll;i<curr_exchange.m_all_num_tiledb_row_idx_vec_request[partition_idx];++i)
  {
    assert(static_cast<size_t>(begin_idx+i) < all_partitions_tiledb_row_idx_vec.size());
    auto row_idx = all_partitions_tiledb_row_idx_vec[begin_idx+i];
    if(row_idx < 0)
      break; 
    int64_t local_file_idx = -1;
    m_vid_mapper->get_local_file_idx_for_row(row_idx, local_file_idx);
    assert(local_file_idx >= 0 && static_cast<size_t>(local_file_idx) < m_vcf2binary_handlers.size());
    m_partition_batch[partition_idx].activate_file(local_file_idx);
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
  for(auto& vcf_handler : m_vcf2binary_handlers)
    vcf_handler.read_next_batch(m_cell_data_buffers, m_partition_batch, false);
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
        int64_t local_file_idx = -1;
        auto status = 
          m_vid_mapper->get_local_file_idx_for_row(curr_exchange.m_all_tiledb_row_idx_vec_response[idx_offset+i],
              local_file_idx);
        assert(status);
        auto& file_batch = m_partition_batch[partition_idx].get_partition_file_batch(local_file_idx);
        auto offset = m_partition_batch[partition_idx].get_partition_begin_offset() + i*m_max_size_per_callset;
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
  m_vid_mapper = static_cast<VidMapper*>(new FileBasedVidMapper(m_vid_mapping_filename));
  m_max_size_per_callset = m_per_partition_size/m_vid_mapper->get_num_callsets();
  //Allocate buffers
  for(auto i=0u;i<m_ping_pong_buffers.size();++i)
    m_ping_pong_buffers[i].resize(m_per_partition_size);
  //Converter processes run independent of loader when num_converter_processes > 0
  if(m_standalone_converter_process)
  {
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
}

void VCF2TileDBLoader::read_all()
{
  while(true)
  {
    m_converter->read_next_batch(0u);
    if(m_owned_exchanges[0].is_new_data_in_converter_response(0u))
      m_converter->dump_latest_buffer(0u, std::cout);
    else
      break;
  }
}

void VCF2TileDBLoader::clear()
{
}
