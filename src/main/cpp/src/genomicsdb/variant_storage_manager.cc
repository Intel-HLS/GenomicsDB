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

#include "variant_storage_manager.h"
#include "variant_field_data.h"
#include <sys/stat.h>
#include "json_config.h"
#include <dirent.h>
#include <uuid/uuid.h>

#define VERIFY_OR_THROW(X) if(!(X)) throw VariantStorageManagerException(#X);
#define METADATA_FILE_PREFIX "genomicsdb_meta"
#define GET_OLD_METADATA_PATH(workspace, array) ((workspace)+'/'+(array)+'/' + METADATA_FILE_PREFIX + ".json")
#define GET_METADATA_DIRECTORY(workspace, array) ((workspace)+'/'+(array)+"/genomicsdb_meta_dir")
inline bool is_genomicsdb_meta_file(const std::string& filename, const std::string& filepath,
    const unsigned char filetype)
{
  return (!filename.empty() && !filepath.empty() && filetype == DT_REG
      && filename.find(METADATA_FILE_PREFIX) != std::string::npos && filename.find(".json") != std::string::npos);
}

const std::unordered_map<std::string, int> VariantStorageManager::m_mode_string_to_int = {
  { "r", TILEDB_ARRAY_READ },
  { "w", TILEDB_ARRAY_WRITE },
  { "a", TILEDB_ARRAY_WRITE }
};

std::vector<const char*> VariantStorageManager::m_metadata_attributes = std::vector<const char*>({ "max_valid_row_idx_in_array" });

//ceil(buffer_size/field_size)*field_size
#define GET_ALIGNED_BUFFER_SIZE(buffer_size, field_size) ((((buffer_size)+(field_size)-1u)/(field_size))*(field_size))

//VariantArrayCellIterator functions
VariantArrayCellIterator::VariantArrayCellIterator(TileDB_CTX* tiledb_ctx, const VariantArraySchema& variant_array_schema,
        const std::string& array_path, const int64_t* range, const std::vector<int>& attribute_ids, const size_t buffer_size)
  : m_num_queried_attributes(attribute_ids.size()), m_tiledb_ctx(tiledb_ctx),
  m_variant_array_schema(&variant_array_schema), m_cell(variant_array_schema, attribute_ids)
#ifdef DO_PROFILING
  , m_tiledb_timer()
  , m_tiledb_to_buffer_cell_timer()
#endif
{
  m_buffers.clear();
  std::vector<const char*> attribute_names(attribute_ids.size()+1u);  //+1 for the COORDS
  for(auto i=0ull;i<attribute_ids.size();++i)
  {
    //Buffer size must be resized to be a multiple of the field size
    auto curr_buffer_size = buffer_size;
    attribute_names[i] = variant_array_schema.attribute_name(attribute_ids[i]).c_str();
    //For varible length attributes, need extra buffer for maintaining offsets
    if(variant_array_schema.is_variable_length_field(attribute_ids[i]))
    {
      curr_buffer_size = GET_ALIGNED_BUFFER_SIZE(buffer_size, sizeof(size_t));
      m_buffers.emplace_back(curr_buffer_size);
    }
    else
      curr_buffer_size = GET_ALIGNED_BUFFER_SIZE(buffer_size, m_cell.get_field_size_in_bytes(i));
    m_buffers.emplace_back(curr_buffer_size);
  }
  //Co-ordinates
  attribute_names[attribute_ids.size()] = TILEDB_COORDS;
  m_buffers.emplace_back(GET_ALIGNED_BUFFER_SIZE(buffer_size, variant_array_schema.dim_size_in_bytes()));
  //Initialize pointers to buffers
  m_buffer_pointers.resize(m_buffers.size());
  m_buffer_sizes.resize(m_buffers.size());
  for(auto i=0ull;i<m_buffers.size();++i)
  {
    m_buffer_pointers[i] = reinterpret_cast<void*>(&(m_buffers[i][0]));
    m_buffer_sizes[i] = m_buffers[i].size();
  }
  /* Initialize the array in READ mode. */
  auto status = tiledb_array_iterator_init(
      tiledb_ctx, 
      &m_tiledb_array_iterator,
      array_path.c_str(),
      TILEDB_ARRAY_READ,
      reinterpret_cast<const void*>(range), // range, 
      &(attribute_names[0]),           
      attribute_names.size(),
      const_cast<void**>(&(m_buffer_pointers[0])),
      &(m_buffer_sizes[0]));
  if(status != TILEDB_OK)
    throw VariantStorageManagerException(std::string("Error while initializing TileDB iterator")
          +"\nTileDB error message : "+tiledb_errmsg);
#ifdef DEBUG
  m_last_row = -1;
  m_last_column = -1;
  m_num_cells_iterated_over = 0ull;
#endif
#ifdef DO_PROFILING
  m_tiledb_timer.stop();
#endif
}

const BufferVariantCell& VariantArrayCellIterator::operator*()
{
#ifdef DO_PROFILING
  m_tiledb_to_buffer_cell_timer.start();
#endif
  const uint8_t* field_ptr = 0;
  size_t field_size = 0u;
  for(auto i=0u;i<m_num_queried_attributes;++i)
  {
    auto status = tiledb_array_iterator_get_value(m_tiledb_array_iterator, i,
        reinterpret_cast<const void**>(&field_ptr), &field_size);
    if(status != TILEDB_OK)
      throw VariantStorageManagerException(std::string("Error while getting value from TileDB iterator")
          +" for field with query index "+std::to_string(i)
          +"\nTileDB error message : "+tiledb_errmsg);
    m_cell.set_field_ptr_for_query_idx(i, field_ptr);
    m_cell.set_field_size_in_bytes(i, field_size);
  }
  //Co-ordinates
  auto status = tiledb_array_iterator_get_value(m_tiledb_array_iterator, m_num_queried_attributes,
        reinterpret_cast<const void**>(&field_ptr), &field_size);
  if(status != TILEDB_OK)
    throw VariantStorageManagerException(std::string("Error while getting coordinate value from TileDB iterator")
        +"\nTileDB error message : "+tiledb_errmsg);
  assert(field_size == m_variant_array_schema->dim_size_in_bytes());
  auto coords_ptr = reinterpret_cast<const int64_t*>(field_ptr);
  m_cell.set_coordinates(coords_ptr[0], coords_ptr[1]);
#ifdef DO_PROFILING
  m_tiledb_to_buffer_cell_timer.stop();
#endif
  return m_cell;
}

//VariantArrayInfo functions
VariantArrayInfo::VariantArrayInfo(int idx, int mode,
    const std::string& workspace, const std::string& name,
    const VidMapper* vid_mapper,
    const VariantArraySchema& schema,
    TileDB_CTX* tiledb_ctx,
    TileDB_Array* tiledb_array,
    const size_t buffer_size)
: m_idx(idx), m_mode(mode),
  m_workspace(workspace), m_name(name),
  m_vid_mapper(vid_mapper), m_schema(schema),m_tiledb_ctx(tiledb_ctx), m_cell(m_schema), m_tiledb_array(tiledb_array)
{
  //If writing, allocate buffers
  if(mode == TILEDB_ARRAY_WRITE || mode == TILEDB_ARRAY_WRITE_UNSORTED)
  {
    m_buffers.clear();
    for(auto i=0ull;i<schema.attribute_num();++i)
    {
      //For varible length attributes, need extra buffer for maintaining offsets
      if(m_schema.is_variable_length_field(i))
        m_buffers.emplace_back(buffer_size);
      m_buffers.emplace_back(buffer_size);
    }
    //Co-ordinates
    m_buffers.emplace_back(buffer_size);
    //Initialize pointers to buffers
    m_buffer_pointers.resize(m_buffers.size());
    m_buffer_offsets.resize(m_buffers.size());
    for(auto i=0ull;i<m_buffers.size();++i)
    {
      m_buffer_pointers[i] = reinterpret_cast<void*>(&(m_buffers[i][0]));
      m_buffer_offsets[i] = 0ull; //will be modified during a write
    }
  }
  read_row_bounds_from_metadata();
#ifdef DEBUG
  m_last_row = m_last_column = -1;
#endif
}

//Move constructor
VariantArrayInfo::VariantArrayInfo(VariantArrayInfo&& other)
  : m_schema(std::move(other.m_schema)), m_cell(std::move(other.m_cell))
{
  m_idx = other.m_idx;
  m_mode = other.m_mode;
  m_workspace = std::move(other.m_workspace);
  m_name = std::move(other.m_name);
  //Pointer handling
  m_tiledb_array = other.m_tiledb_array;
  other.m_tiledb_array = 0;
  m_vid_mapper = other.m_vid_mapper;
  other.m_vid_mapper = 0;
  //Point array schema to this array schema
  m_cell.set_variant_array_schema(m_schema);
  //Move other members
  m_buffers = std::move(other.m_buffers);
  m_buffer_offsets = std::move(other.m_buffer_offsets);
  m_buffer_pointers = std::move(other.m_buffer_pointers);
  for(auto i=0ull;i<m_buffer_pointers.size();++i)
    m_buffer_pointers[i] = reinterpret_cast<void*>(&(m_buffers[i][0]));
  m_metadata_contains_max_valid_row_idx_in_array = other.m_metadata_contains_max_valid_row_idx_in_array;
  m_max_valid_row_idx_in_array = other.m_max_valid_row_idx_in_array;
  m_metadata_contains_lb_row_idx = other.m_metadata_contains_lb_row_idx;
  m_lb_row_idx = other.m_lb_row_idx;
#ifdef DEBUG
  m_last_row = other.m_last_row;
  m_last_column = other.m_last_column;
#endif
}

void VariantArrayInfo::write_cell(const void* ptr)
{
  m_cell.set_cell(ptr);
#ifdef DEBUG
  assert((m_cell.get_begin_column() > m_last_column) || (m_cell.get_begin_column() == m_last_column && m_cell.get_row() > m_last_row));
  m_last_row = m_cell.get_row();
  m_last_column = m_cell.get_begin_column();
#endif
  auto buffer_idx = 0ull;
  auto overflow = false;
  //First check if the current cell will fit into the buffers
  for(auto i=0ull;i<m_schema.attribute_num();++i)
  {
    assert(buffer_idx < m_buffer_pointers.size());
    //Specify offsets in the buffer for variable length fields
    if(m_schema.is_variable_length_field(i))
    {
      if(m_buffer_offsets[buffer_idx]+sizeof(size_t) > m_buffers[buffer_idx].size())
      {
        overflow = true;
        break;
      }
      ++buffer_idx;
    }
    if(m_buffer_offsets[buffer_idx]+m_cell.get_field_size_in_bytes(i) > m_buffers[buffer_idx].size())
    {
      overflow = true;
      break;
    }
    ++buffer_idx;
  }
  //Check if co-ordinates buffer overflows
  auto coords_buffer_idx = m_buffers.size()-1u;
  auto coords_size = m_schema.dim_size_in_bytes();
  overflow = overflow || (m_buffer_offsets[coords_buffer_idx]+coords_size > m_buffers[coords_buffer_idx].size());
  //write to array and reset sizes
  if(overflow)
  {
    auto status = tiledb_array_write(m_tiledb_array, const_cast<const void**>(&(m_buffer_pointers[0])), &(m_buffer_offsets[0]));
    if(status != TILEDB_OK)
      throw VariantStorageManagerException(std::string("Error while writing to TileDB array")
          +"\nTileDB error message : "+tiledb_errmsg);
    memset(&(m_buffer_offsets[0]), 0, m_buffer_offsets.size()*sizeof(size_t));
  }
  buffer_idx = 0;
  for(auto i=0ull;i<m_schema.attribute_num();++i)
  {
    assert(buffer_idx < m_buffer_pointers.size());
    //Specify offsets in the buffer for variable length fields
    if(m_schema.is_variable_length_field(i))
    {
      assert(buffer_idx+1u < m_buffer_offsets.size());
      //This could happen if segment_size in the loader JSON is so small that the data for a single cell
      //does not fit into the buffer
      if(overflow && m_buffer_offsets[buffer_idx]+sizeof(size_t) > m_buffers[buffer_idx].size())
      {
        m_buffers[buffer_idx].resize(m_buffer_offsets[buffer_idx]+sizeof(size_t));
        m_buffer_pointers[buffer_idx] = reinterpret_cast<void*>(&(m_buffers[buffer_idx][0u]));
      }
      assert(m_buffer_offsets[buffer_idx]+sizeof(size_t) <= m_buffers[buffer_idx].size());
      //Offset buffer - new entry starts after last entry
      *(reinterpret_cast<size_t*>(&(m_buffers[buffer_idx][m_buffer_offsets[buffer_idx]]))) = m_buffer_offsets[buffer_idx+1u];
      m_buffer_offsets[buffer_idx] += sizeof(size_t);
      ++buffer_idx;
    }
    auto field_size = m_cell.get_field_size_in_bytes(i);
    //This could happen if segment_size in the loader JSON is so small that the data for a single cell
    //does not fit into the buffer
    if(overflow && m_buffer_offsets[buffer_idx]+field_size > m_buffers[buffer_idx].size())
    {
      m_buffers[buffer_idx].resize(m_buffer_offsets[buffer_idx]+field_size);
      m_buffer_pointers[buffer_idx] = reinterpret_cast<void*>(&(m_buffers[buffer_idx][0u]));
    }
    assert(m_buffer_offsets[buffer_idx]+field_size <= m_buffers[buffer_idx].size());
    memcpy(&(m_buffers[buffer_idx][m_buffer_offsets[buffer_idx]]), m_cell.get_field_ptr_for_query_idx<void>(i), field_size);
    m_buffer_offsets[buffer_idx] += field_size;
    ++buffer_idx;
  }
  //Co-ordinates
  assert(buffer_idx == coords_buffer_idx);
  assert(m_buffer_offsets[coords_buffer_idx]+coords_size <= m_buffers[coords_buffer_idx].size());
  memcpy(&(m_buffers[coords_buffer_idx][m_buffer_offsets[coords_buffer_idx]]), ptr, coords_size);
  m_buffer_offsets[coords_buffer_idx] += coords_size;
}

void VariantArrayInfo::read_row_bounds_from_metadata()
{
  //Compute value from array schema
  m_metadata_contains_max_valid_row_idx_in_array = false;
  m_metadata_contains_lb_row_idx = false;
  const auto& dim_domains = m_schema.dim_domains();
  m_lb_row_idx = dim_domains[0].first;
  m_max_valid_row_idx_in_array = dim_domains[0].second;
  //To read genomicsdb json entries in meta directory
  std::string buffer;
  std::vector<char*> dir_entries_pointers;
  std::vector<unsigned char> dir_entries_type;
  get_genomicsdb_meta_directory_entries(buffer, dir_entries_pointers, dir_entries_type);
  auto metadata_directory = GET_METADATA_DIRECTORY(m_workspace, m_name);
  for(auto i=0u;i<dir_entries_pointers.size();++i)
    read_row_bounds_from_metadata(dir_entries_pointers[i],
        metadata_directory + '/' + dir_entries_pointers[i],
        dir_entries_type[i]);
  //For the old metadata file
  read_row_bounds_from_metadata(std::string(METADATA_FILE_PREFIX)+".json",
      GET_OLD_METADATA_PATH(m_workspace, m_name),
      DT_REG);
}

int VariantArrayInfo::read_row_bounds_from_metadata(const std::string& filename, const std::string& filepath,
    const unsigned char filetype)
{
  //Try reading from metadata
  if(is_genomicsdb_meta_file(filename, filepath, filetype))
  {
    std::ifstream ifs(filepath.c_str());
    if(ifs.is_open())
    {
      rapidjson::Document json_doc;
      std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
      json_doc.Parse(str.c_str());
      if(json_doc.HasParseError())
        return -1;
        //throw VariantStorageManagerException(std::string("Syntax error in corrupted JSON metadata file ")+m_metadata_filename);
      if(json_doc.HasMember("max_valid_row_idx_in_array") && json_doc["max_valid_row_idx_in_array"].IsInt64()
          && (!m_metadata_contains_max_valid_row_idx_in_array ||
           m_max_valid_row_idx_in_array < json_doc["max_valid_row_idx_in_array"].GetInt64())
          )
      {
        m_max_valid_row_idx_in_array = json_doc["max_valid_row_idx_in_array"].GetInt64();
        m_metadata_contains_max_valid_row_idx_in_array = true;
      }
      if(json_doc.HasMember("lb_row_idx") && json_doc["lb_row_idx"].IsInt64()
          && (!m_metadata_contains_lb_row_idx ||
            m_lb_row_idx > json_doc["lb_row_idx"].GetInt64()))
      {
        m_lb_row_idx = json_doc["lb_row_idx"].GetInt64();
        m_metadata_contains_lb_row_idx = true;
      }
      return 0;
    }
  }
  return -1;
}

void VariantArrayInfo::update_row_bounds_in_array(
    const int64_t lb_row_idx, const int64_t max_valid_row_idx_in_array)
{
  auto update_metadata = false;
  //Update metadata if:

  // (it did not exist and (#valid rows is set to large value or  num_rows_seen > num rows as defined in schema
  //                (old implementation))  OR
  // (num_rows_seen > #valid rows in metadata (this part implies that metadata file exists)
  if((!m_metadata_contains_max_valid_row_idx_in_array &&
        (m_max_valid_row_idx_in_array == INT64_MAX-1 || max_valid_row_idx_in_array > m_max_valid_row_idx_in_array))
      || (max_valid_row_idx_in_array > m_max_valid_row_idx_in_array))
  {
    m_max_valid_row_idx_in_array = max_valid_row_idx_in_array;
    update_metadata = true;
  }
  // (metadata did not exist OR
  // (new_lb < value in metadata file)
  if(!m_metadata_contains_lb_row_idx || lb_row_idx < m_lb_row_idx)
  {
    m_lb_row_idx = lb_row_idx;
    update_metadata = true;
  }
  if(update_metadata)
  {
    rapidjson::Document json_doc;
    json_doc.SetObject();
    json_doc.AddMember("lb_row_idx", m_lb_row_idx, json_doc.GetAllocator());
    json_doc.AddMember("max_valid_row_idx_in_array", m_max_valid_row_idx_in_array, json_doc.GetAllocator());
    rapidjson::StringBuffer buffer;
    rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
    json_doc.Accept(writer);
    uuid_t value;
    uuid_generate(value);
    char uuid_str[40];
    uuid_unparse(value, uuid_str);
    auto metadata_filepath = GET_METADATA_DIRECTORY(m_workspace, m_name)
      + '/' + METADATA_FILE_PREFIX
      + '_' + uuid_str +  ".json";
    auto fptr = tiledb_open_file(metadata_filepath.c_str(), "w");
    if(fptr == 0)
      throw VariantStorageManagerException(std::string("Could not open metadata file ")+metadata_filepath
          +"\nTileDB error message : "+tiledb_errmsg);
    auto string_length = strlen(buffer.GetString());
    auto num_bytes_written =
      tiledb_fwrite(reinterpret_cast<const void*>(buffer.GetString()), 1u, string_length, fptr);
    if(num_bytes_written != string_length)
      throw VariantStorageManagerException(std::string("Could not write to metadata file ")+metadata_filepath
          +"\nTileDB error message : "+tiledb_errmsg);
    tiledb_close_file(fptr);
  }
}

void VariantArrayInfo::close_array(const bool consolidate_tiledb_array)
{
  auto coords_buffer_idx = m_buffers.size()-1u;
  if((m_mode == TILEDB_ARRAY_WRITE || m_mode == TILEDB_ARRAY_WRITE_UNSORTED))
  {
    //Flush cells in buffer
    if(m_buffer_offsets[coords_buffer_idx] > 0ull)
    {
      auto status = tiledb_array_write(m_tiledb_array, const_cast<const void**>(&(m_buffer_pointers[0])), &(m_buffer_offsets[0]));
      if(status != TILEDB_OK)
        throw VariantStorageManagerException("Error while writing to array "+m_name
            +"\nTileDB error message : "+tiledb_errmsg);
      memset(&(m_buffer_offsets[0]), 0, m_buffer_offsets.size()*sizeof(size_t));
    }
    auto sync_status = tiledb_array_sync(m_tiledb_array);
    if(sync_status != TILEDB_OK)
      throw VariantStorageManagerException("Error while syncing array "+m_name+" to disk"
          +"\nTileDB error message : "+tiledb_errmsg);
  }
  if(m_tiledb_array)
  {
    auto status = tiledb_array_finalize(m_tiledb_array);
    if(status != TILEDB_OK)
      throw VariantStorageManagerException("Error while finalizing TileDB array "+m_name
          +"\nTileDB error message : "+tiledb_errmsg);
    if(consolidate_tiledb_array)
    {
      //Consolidate metadb entries as well
      std::string buffer;
      std::vector<char*> dir_entries_pointers;
      std::vector<unsigned char> dir_entries_type;
      get_genomicsdb_meta_directory_entries(buffer, dir_entries_pointers, dir_entries_type);
      auto metadata_directory = GET_METADATA_DIRECTORY(m_workspace, m_name);
      for(auto i=0u;i<dir_entries_pointers.size();++i)
      {
        std::string filename = dir_entries_pointers[i];
        std::string filepath = metadata_directory+'/'+dir_entries_pointers[i];
        if(is_genomicsdb_meta_file(filename, filepath, dir_entries_type[i]))
        {
          auto read_status = read_row_bounds_from_metadata(filename, filepath, dir_entries_type[i]);
          if(read_status == 0) //no errors
            tiledb_delete_file(filepath.c_str());
        }
      }
      //Force write out consolidated values
      auto lb_row_idx = m_lb_row_idx;
      m_lb_row_idx = INT64_MAX;
      auto max_valid_row_idx_in_array = m_max_valid_row_idx_in_array;
      m_max_valid_row_idx_in_array = -1;
      update_row_bounds_in_array(lb_row_idx, max_valid_row_idx_in_array);
      auto status = tiledb_array_consolidate(m_tiledb_ctx, (m_workspace + '/' + m_name).c_str());
      if(status != TILEDB_OK)
        throw VariantStorageManagerException("Error while consolidating TileDB array "+m_name
            +"\nTileDB error message : "+tiledb_errmsg);
    }
  }
  m_tiledb_array = 0;
  m_name.clear();
  m_mode = -1;
}

void VariantArrayInfo::get_genomicsdb_meta_directory_entries(std::string& buffer,
    std::vector<char*>& dir_entries_pointers,
    std::vector<unsigned char>& dir_entries_type)
{
  auto metadata_directory = GET_METADATA_DIRECTORY(m_workspace, m_name);
  //To read genomicsdb json entries in meta directory
  size_t num_entries_in_directory = 4096u;
  auto retry = true;
  while(retry)
  {
    buffer.resize((NAME_MAX+1u)*num_entries_in_directory);
    dir_entries_pointers.resize(num_entries_in_directory);
    dir_entries_type.resize(num_entries_in_directory);
    for(auto i=0u;i<num_entries_in_directory;++i)
      dir_entries_pointers[i] = &(buffer[i*(NAME_MAX+1u)]);
    tiledb_ls_directory(metadata_directory.c_str(),
        &(dir_entries_pointers[0u]),
        &(dir_entries_type[0u]),
        &num_entries_in_directory);
    if(num_entries_in_directory < dir_entries_pointers.size()) //no overflow
      retry = false;
    else
    {
      if(num_entries_in_directory > 100000000u) //100M
        throw VariantStorageManagerException(std::string("GenomicsDB meta directory ")
            +metadata_directory+" has too many entries - cannot handle so many entries");
      num_entries_in_directory *= 2u; //retry with double the number of entries
    }
  }
  dir_entries_pointers.resize(num_entries_in_directory);
  dir_entries_type.resize(num_entries_in_directory);
}

//VariantStorageManager functions
VariantStorageManager::VariantStorageManager(const std::string& workspace, const unsigned segment_size)
{
  m_workspace = workspace;
  m_segment_size = segment_size;
  /*Initialize context with default params*/
  tiledb_ctx_init(&m_tiledb_ctx, NULL);
  //Create workspace if it does not exist
  struct stat st;
  auto status = stat(workspace.c_str(), &st);
  //Exists and is not a directory
  if(status >= 0 && !S_ISDIR(st.st_mode))
    throw VariantStorageManagerException(std::string("Workspace path ")+workspace+" exists and is not a directory");
  //Doesn't exist, create workspace
  if(status < 0)
    if(tiledb_workspace_create(m_tiledb_ctx, workspace.c_str()) != TILEDB_OK)
      throw VariantStorageManagerException(std::string("Error while creating TileDB workspace ")
          + workspace
          +"\nTileDB error message : "+tiledb_errmsg);
}

bool VariantStorageManager::check_if_TileDB_array_exists(const std::string& array_name)
{
  //Use tiledb_ls call to avoid non-faulty error message while trying to init array during loading
  std::vector<char*> ws_entries(1024u);
  std::vector<int> ws_entry_types;
  for(auto i=0ull;i<ws_entries.size();++i)
    ws_entries[i] = new char[TILEDB_NAME_MAX_LEN+1]; //for null char
  auto num_entries = 0;
  auto ls_status = TILEDB_ERR;
  while(ls_status == TILEDB_ERR && ws_entries.size() < 100000ll) //cutoff if too many entries in workspace
  {
    num_entries = ws_entries.size();
    ws_entry_types.resize(ws_entries.size());
    ls_status = tiledb_ls(m_tiledb_ctx, m_workspace.c_str(), &(ws_entries[0]), &(ws_entry_types[0]), &num_entries);
    if(ls_status == TILEDB_ERR)
    {
      std::cerr << "[GenomicsDB::VariantStorageManager] INFO: ignore message \"[TileDB::StorageManager] Error: Cannot list TileDB directory; Directory buffer overflow.\" in the previous line\n";
      auto old_size = ws_entries.size();
      ws_entries.resize(2u*ws_entries.size()+1u);
      for(auto i=old_size;i<ws_entries.size();++i)
        ws_entries[i] = new char[TILEDB_NAME_MAX_LEN+1]; //for null char
    }
  }
  if(ls_status == TILEDB_ERR)
    throw VariantStorageManagerException(std::string("Too many entries in the workspace ")+m_workspace+" - cannot handle more than 100K entries");
  auto string_length = std::min<size_t>(array_name.length(), TILEDB_NAME_MAX_LEN);
  auto array_exists = false;
  for(auto i=0ull;i<static_cast<size_t>(num_entries);++i)
  {
    if(ws_entry_types[i] == TILEDB_ARRAY
        && strnlen(ws_entries[i], TILEDB_NAME_MAX_LEN+1u) == string_length
        && strncmp(array_name.c_str(), ws_entries[i], string_length) == 0)
      array_exists = true;
    delete[] ws_entries[i];
  }
  for(auto i=static_cast<size_t>(num_entries);i<ws_entries.size();++i)
    delete[] ws_entries[i];
  return array_exists;
}

int VariantStorageManager::open_array(const std::string& array_name, const VidMapper* vid_mapper, const char* mode,
    const int tiledb_zlib_compression_level)
{
  auto mode_iter = VariantStorageManager::m_mode_string_to_int.find(mode);
  VERIFY_OR_THROW(mode_iter != VariantStorageManager::m_mode_string_to_int.end() && "Unknown mode of opening an array");
  auto mode_int = (*mode_iter).second;
  auto array_exists = check_if_TileDB_array_exists(array_name);
  if(array_exists)
  {
    //Try to open the array
    TileDB_Array* tiledb_array;
    auto status = tiledb_array_init(
        m_tiledb_ctx, 
        &tiledb_array,
        (m_workspace+'/'+array_name).c_str(),
        mode_int, NULL,
        0, 0);
    if(status == TILEDB_OK)
    {
      auto idx = m_open_arrays_info_vector.size();
      //Schema
      VariantArraySchema tmp_schema;
      get_array_schema(array_name, &tmp_schema);
      //Just to be safe
      define_metadata_schema(&tmp_schema);
      m_open_arrays_info_vector.emplace_back(idx, mode_int,
          m_workspace, array_name,
          vid_mapper, tmp_schema,
          m_tiledb_ctx, tiledb_array,
          m_segment_size);
      tiledb_array_set_zlib_compression_level(tiledb_array, tiledb_zlib_compression_level);
      return idx;
    }
  }
  return -1;
}

void VariantStorageManager::close_array(const int ad, const bool consolidate_tiledb_array)
{
  VERIFY_OR_THROW(static_cast<size_t>(ad) < m_open_arrays_info_vector.size() &&
      m_open_arrays_info_vector[ad].get_array_name().length());
  m_open_arrays_info_vector[ad].close_array(consolidate_tiledb_array);
}

int VariantStorageManager::define_array(const VariantArraySchema* variant_array_schema, const size_t num_cells_per_tile)
{
  //Attribute attributes
  std::vector<const char*> attribute_names(variant_array_schema->attribute_num());
  std::vector<int> cell_val_num(variant_array_schema->attribute_num());
  std::vector<int> types(variant_array_schema->attribute_num()+1u);  //+1 for the co-ordinates
  std::vector<int> compression(variant_array_schema->attribute_num()+1u);  //+1 for the co-ordinates
  for(auto i=0ull;i<variant_array_schema->attribute_num();++i)
  {
    attribute_names[i] = variant_array_schema->attribute_name(i).c_str();
    cell_val_num[i] = variant_array_schema->val_num(i);
    types[i] = VariantFieldTypeUtil::get_tiledb_type_for_variant_field_type(variant_array_schema->type(i));
    compression[i] = variant_array_schema->compression(i);
  }
  //Co-ordinates
  types[variant_array_schema->attribute_num()] = VariantFieldTypeUtil::get_tiledb_type_for_variant_field_type(
      variant_array_schema->dim_type());
  compression[variant_array_schema->attribute_num()] = variant_array_schema->dim_compression_type();
  std::vector<const char*> dim_names(variant_array_schema->dim_names().size());
  std::vector<int64_t> dim_domains(2u*dim_names.size());
  for(auto i=0ull;i<dim_names.size();++i)
  {
    dim_names[i] = variant_array_schema->dim_names()[i].c_str();
    dim_domains[2u*i] = variant_array_schema->dim_domains()[i].first;
    dim_domains[2u*i+1u] = variant_array_schema->dim_domains()[i].second;
  }
  //TileDB C API
  TileDB_ArraySchema array_schema;
  memset(&array_schema, 0, sizeof(TileDB_ArraySchema));
  tiledb_array_set_schema(
      // The array schema struct
      &array_schema,
      // Array name
      (m_workspace+'/'+variant_array_schema->array_name()).c_str(),
      // Attributes
      &(attribute_names[0]),
      // Number of attributes
      attribute_names.size(),
      // Capacity
      num_cells_per_tile,
      // Cell order
      TILEDB_COL_MAJOR,
      // Number of cell values per attribute (NULL means 1 everywhere)
      &(cell_val_num[0]),
      // Compression
      &(compression[0]),
      // Sparse array
      0,
      // Dimensions
      &(dim_names[0]),
      // Number of dimensions
      dim_names.size(),
       // Domain
      &(dim_domains[0]),
      // Domain length in bytes
      dim_domains.size()*sizeof(int64_t),
      // Tile extents (no regular tiles defined)
      NULL,
      // Tile extents in bytes
      0, 
      // Tile order (0 means ignore in sparse arrays and default in dense)
      0,
      // Types
      &(types[0])
  );
  /* Create the array schema */
  auto status = tiledb_array_create(m_tiledb_ctx, &array_schema);
  if(status == TILEDB_OK)
  {
    status = tiledb_array_free_schema(&array_schema);
    if(status == TILEDB_OK)
      status = define_metadata_schema(variant_array_schema);
  }
  return status;
}

void VariantStorageManager::delete_array(const std::string& array_name)
{
  auto array_exists = check_if_TileDB_array_exists(array_name);
  if(array_exists)
  {
    tiledb_delete_file(GET_OLD_METADATA_PATH(m_workspace, array_name).c_str());
    tiledb_delete_directory(GET_METADATA_DIRECTORY(m_workspace, array_name).c_str());
    auto status = tiledb_delete(m_tiledb_ctx, (m_workspace+"/"+array_name).c_str());
    if(status != TILEDB_OK)
      throw VariantStorageManagerException(std::string("Error while deleting TileDB array ")
          + (m_workspace+"/"+array_name)
          + "\nTileDB error message : "+tiledb_errmsg);
  }
}

//Define metadata
int VariantStorageManager::define_metadata_schema(const VariantArraySchema* variant_array_schema)
{
  auto status = tiledb_create_directory(GET_METADATA_DIRECTORY(m_workspace, variant_array_schema->array_name()).c_str());
  if(status != TILEDB_OK)
    throw VariantStorageManagerException(std::string("Could not create GenomicsDB meta-data directory ")
          + "\nTileDB error message : "+tiledb_errmsg);
  return TILEDB_OK;
}

int VariantStorageManager::get_array_schema(const int ad, VariantArraySchema* variant_array_schema)
{
  VERIFY_OR_THROW(static_cast<size_t>(ad) < m_open_arrays_info_vector.size() &&
      m_open_arrays_info_vector[ad].get_array_name().length());
  auto status = get_array_schema(m_open_arrays_info_vector[ad].get_array_name(), variant_array_schema);
  if(status == TILEDB_OK)
    m_open_arrays_info_vector[ad].set_schema(*variant_array_schema);
  return status;
}

int VariantStorageManager::get_array_schema(const std::string& array_name, VariantArraySchema* variant_array_schema)
{
  // Get array schema
  TileDB_ArraySchema tiledb_array_schema;
  memset(&tiledb_array_schema, 0, sizeof(TileDB_ArraySchema));
  auto status = tiledb_array_load_schema(m_tiledb_ctx, (m_workspace+'/'+array_name).c_str(),
    &tiledb_array_schema);
  if(status != TILEDB_OK)
    return -1;
  //Attribute attributes
  std::vector<std::string> attribute_names(tiledb_array_schema.attribute_num_);
  std::vector<int> val_num(tiledb_array_schema.attribute_num_);
  std::vector<std::type_index> attribute_types(tiledb_array_schema.attribute_num_+1u, std::type_index(typeid(void)));//+1 for co-ordinates
  std::vector<int> compression(tiledb_array_schema.attribute_num_+1u);  //+1 for co-ordinates
  for(auto i=0u;i<attribute_names.size();++i)
  {
    attribute_names[i] = tiledb_array_schema.attributes_[i];
    val_num[i] = tiledb_array_schema.cell_val_num_[i];
    attribute_types[i] = VariantFieldTypeUtil::get_variant_field_type_for_tiledb_type(tiledb_array_schema.types_[i]);
    compression[i] = tiledb_array_schema.compression_[i];
  }
  //Co-ords
  auto coords_idx = tiledb_array_schema.attribute_num_;
  attribute_types[coords_idx] = std::type_index(typeid(int64_t));
  compression[coords_idx] = tiledb_array_schema.compression_[coords_idx];
  std::vector<std::string> dim_names(tiledb_array_schema.dim_num_);
  auto dim_domains = std::vector<std::pair<int64_t,int64_t>>(tiledb_array_schema.dim_num_);
  auto dim_domains_ptr = reinterpret_cast<const int64_t*>(tiledb_array_schema.domain_);
  for(auto i=0;i<tiledb_array_schema.dim_num_;++i)
  {
    dim_names[i] = tiledb_array_schema.dimensions_[i];
    dim_domains[i].first = dim_domains_ptr[2*i];
    dim_domains[i].second = dim_domains_ptr[2*i+1];
  }
  *variant_array_schema = std::move(VariantArraySchema(
        array_name,
        attribute_names,
        dim_names,
        dim_domains,
        attribute_types,
        val_num, 
        compression,
        TILEDB_COL_MAJOR));
  // Free array schema
  tiledb_array_free_schema(&tiledb_array_schema);
  return TILEDB_OK;
}

VariantArrayCellIterator* VariantStorageManager::begin(
    int ad, const int64_t* range, const std::vector<int>& attribute_ids) const
{
  VERIFY_OR_THROW(static_cast<size_t>(ad) < m_open_arrays_info_vector.size() &&
      m_open_arrays_info_vector[ad].get_array_name().length());
  auto& curr_elem = m_open_arrays_info_vector[ad];
  return new VariantArrayCellIterator(m_tiledb_ctx, curr_elem.get_schema(), m_workspace+'/'+curr_elem.get_array_name(),
      range, attribute_ids, m_segment_size);   
}

SingleCellTileDBIterator* VariantStorageManager::begin_columnar_iterator(
    int ad, const VariantQueryConfig& query_config, const bool use_common_array_object) const
{
  VERIFY_OR_THROW(static_cast<size_t>(ad) < m_open_arrays_info_vector.size() &&
      m_open_arrays_info_vector[ad].get_array_name().length());
  auto& curr_elem = m_open_arrays_info_vector[ad];
  return new SingleCellTileDBIterator(m_tiledb_ctx,
      use_common_array_object ? m_open_arrays_info_vector[ad].get_tiledb_array() : 0,
      curr_elem.get_vid_mapper(), curr_elem.get_schema(),
      m_workspace+'/'+curr_elem.get_array_name(),
      query_config, m_segment_size);
}

void VariantStorageManager::write_cell_sorted(const int ad, const void* ptr)
{
  assert(static_cast<size_t>(ad) < m_open_arrays_info_vector.size() &&
      m_open_arrays_info_vector[ad].get_array_name().length());
  m_open_arrays_info_vector[ad].write_cell(ptr);
}

int64_t VariantStorageManager::get_num_valid_rows_in_array(const int ad) const
{
  assert(static_cast<size_t>(ad) < m_open_arrays_info_vector.size() &&
      m_open_arrays_info_vector[ad].get_array_name().length());
  return m_open_arrays_info_vector[ad].get_num_valid_rows_in_array();
}

void VariantStorageManager::update_row_bounds_in_array(const int ad, const int64_t lb_row_idx, const int64_t max_valid_row_idx_in_array)
{
  assert(static_cast<size_t>(ad) < m_open_arrays_info_vector.size() &&
      m_open_arrays_info_vector[ad].get_array_name().length());
  m_open_arrays_info_vector[ad].update_row_bounds_in_array(lb_row_idx, max_valid_row_idx_in_array);
}
