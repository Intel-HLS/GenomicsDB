/**
 * Copyright (c) 2016 Intel Corporation
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

#define VERIFY_OR_THROW(X) if(!(X)) throw VariantStorageManagerException(#X);

const std::unordered_map<std::string, int> VariantStorageManager::m_mode_string_to_int = {
  { "r", TILEDB_ARRAY_READ },
  { "w", TILEDB_ARRAY_WRITE },
  { "a", TILEDB_ARRAY_WRITE }
};

//ceil(buffer_size/field_size)*field_size
#define GET_ALIGNED_BUFFER_SIZE(buffer_size, field_size) ((((buffer_size)+(field_size)-1u)/(field_size))*(field_size))

//VariantArrayCellIterator functions
VariantArrayCellIterator::VariantArrayCellIterator(TileDB_CTX* tiledb_ctx, const VariantArraySchema& variant_array_schema,
        const std::string& array_path, const int64_t* range, const std::vector<int>& attribute_ids, const size_t buffer_size)
  : m_num_queried_attributes(attribute_ids.size()), m_tiledb_ctx(tiledb_ctx),
  m_variant_array_schema(&variant_array_schema), m_cell(variant_array_schema, attribute_ids)
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
      reinterpret_cast<const void*>(range), // range, 
      &(attribute_names[0]),           
      attribute_names.size(),
      const_cast<void**>(&(m_buffer_pointers[0])),
      &(m_buffer_sizes[0]));      
  VERIFY_OR_THROW(status == TILEDB_OK && "Error while initializing TileDB iterator");
#ifdef DEBUG
  m_last_row = -1;
  m_last_column = -1;
  m_num_cells_iterated_over = 0ull;
#endif
}

const BufferVariantCell& VariantArrayCellIterator::operator*()
{
  const uint8_t* field_ptr = 0;
  size_t field_size = 0u;
  for(auto i=0u;i<m_num_queried_attributes;++i)
  {
    auto status = tiledb_array_iterator_get_value(m_tiledb_array_iterator, i,
        reinterpret_cast<const void**>(&field_ptr), &field_size);
    VERIFY_OR_THROW(status == TILEDB_OK);
    m_cell.set_field_ptr_for_query_idx(i, field_ptr);
    m_cell.set_field_size_in_bytes(i, field_size);
  }
  //Co-ordinates
  auto status = tiledb_array_iterator_get_value(m_tiledb_array_iterator, m_num_queried_attributes,
        reinterpret_cast<const void**>(&field_ptr), &field_size);
  VERIFY_OR_THROW(status == TILEDB_OK);
  assert(field_size == m_variant_array_schema->dim_size_in_bytes());
  auto coords_ptr = reinterpret_cast<const int64_t*>(field_ptr);
  m_cell.set_coordinates(coords_ptr[0], coords_ptr[1]);
  return m_cell;
}

//VariantArrayInfo functions
VariantArrayInfo::VariantArrayInfo(int idx, int mode, const std::string& name, const VariantArraySchema& schema, TileDB_Array* tiledb_array,
    const size_t buffer_size)
: m_idx(idx), m_mode(mode), m_name(name), m_schema(schema), m_cell(m_schema), m_tiledb_array(tiledb_array)
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
  m_name = std::move(other.m_name);
  //Pointer handling
  m_tiledb_array = other.m_tiledb_array;
  other.m_tiledb_array = 0;
  //Point array schema to this array schema
  m_cell.set_variant_array_schema(m_schema);
  //Move other members
  m_buffers = std::move(other.m_buffers);
  m_buffer_offsets = std::move(other.m_buffer_offsets);
  m_buffer_pointers = std::move(other.m_buffer_pointers);
  for(auto i=0ull;i<m_buffer_pointers.size();++i)
    m_buffer_pointers[i] = reinterpret_cast<void*>(&(m_buffers[i][0]));
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
    VERIFY_OR_THROW(status == TILEDB_OK);
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
      assert(m_buffer_offsets[buffer_idx]+sizeof(size_t) <= m_buffers[buffer_idx].size());
      //Offset buffer - new entry starts after last entry
      *(reinterpret_cast<size_t*>(&(m_buffers[buffer_idx][m_buffer_offsets[buffer_idx]]))) = m_buffer_offsets[buffer_idx+1u];
      m_buffer_offsets[buffer_idx] += sizeof(size_t);
      ++buffer_idx;
    }
    auto field_size = m_cell.get_field_size_in_bytes(i);
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
    VERIFY_OR_THROW(tiledb_workspace_create(m_tiledb_ctx, workspace.c_str()) == TILEDB_OK);
}

int VariantStorageManager::open_array(const std::string& array_name, const char* mode)
{
  auto mode_iter = VariantStorageManager::m_mode_string_to_int.find(mode);
  VERIFY_OR_THROW(mode_iter != VariantStorageManager::m_mode_string_to_int.end() && "Unknown mode of opening an array");
  auto mode_int = (*mode_iter).second;
  //Try to open the array
  TileDB_Array* tiledb_array;
  auto status = tiledb_array_init(
      m_tiledb_ctx, 
      &tiledb_array,
      (m_workspace+'/'+array_name).c_str(),
      mode_int,
      0, 0, 0);
  if(status == TILEDB_OK)
  {
    auto idx = m_open_arrays_info_vector.size();
    VariantArraySchema tmp_schema;
    get_array_schema(array_name, &tmp_schema);
    m_open_arrays_info_vector.emplace_back(idx, mode_int, array_name, tmp_schema, tiledb_array, m_segment_size);
    return idx;
  }
  else
    return -1;
}

void VariantStorageManager::close_array(const int ad)
{
  VERIFY_OR_THROW(static_cast<size_t>(ad) < m_open_arrays_info_vector.size() &&
      m_open_arrays_info_vector[ad].get_array_name().length());
  m_open_arrays_info_vector[ad].close_array();
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
      // Dimensions
      &(dim_names[0]),
      // Number of dimensions
      dim_names.size(),
      // Sparse array
      0,
      // Domain
      &(dim_domains[0]),
      // Domain length in bytes
      dim_domains.size()*sizeof(int64_t),
      // Tile extents (no regular tiles defined)
      NULL,
      // Tile extents in bytes
      0, 
      // Types 
      &(types[0]),
      // Number of cell values per attribute (NULL means 1 everywhere)
      &(cell_val_num[0]),
      // Cell order
      TILEDB_COL_MAJOR,
      // Tile order (0 means ignore in sparse arrays and default in dense)
      0,
      // Capacity
      num_cells_per_tile,
      // Compression
      &(compression[0])
  );
  /* Create the array schema */
  auto status = tiledb_array_create(m_tiledb_ctx, &array_schema);
  if(status == TILEDB_OK)
    status = tiledb_array_free_schema(&array_schema);
  return status;
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

void VariantStorageManager::write_cell_sorted(const int ad, const void* ptr)
{
  assert(static_cast<size_t>(ad) < m_open_arrays_info_vector.size() &&
      m_open_arrays_info_vector[ad].get_array_name().length());
  m_open_arrays_info_vector[ad].write_cell(ptr);
}
