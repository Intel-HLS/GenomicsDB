#include "variant_storage_manager.h"
#include "variant_field_data.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw VariantStorageManagerException(#X);

const std::unordered_map<std::string, int> VariantStorageManager::m_mode_string_to_int = {
  { "r", TILEDB_READ },
  { "w", TILEDB_WRITE }
};

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
    attribute_names[i] = variant_array_schema.attribute_name(attribute_ids[i]).c_str();
    //For varible length attributes, need extra buffer for maintaining offsets
    if(variant_array_schema.is_variable_length_field(attribute_ids[i]))
      m_buffers.emplace_back(buffer_size);
    m_buffers.emplace_back(buffer_size);
  }
  //Co-ordinates
  attribute_names[attribute_ids.size()] = TILEDB_COORDS_NAME;
  m_buffers.emplace_back(buffer_size);
  //Initialize pointers to buffers
  m_buffer_pointers.resize(m_buffers.size());
  m_buffer_sizes.resize(m_buffers.size());
  for(auto i=0ull;i<m_buffers.size();++i)
  {
    m_buffer_pointers[i] = reinterpret_cast<void*>(&(m_buffers[i]));
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
}

const BufferVariantCell& VariantArrayCellIterator::operator*()
{
  const uint8_t* field_ptr = 0;
  size_t field_size = 0u;
  for(auto i=0u;i<m_num_queried_attributes;++i)
  {
    tiledb_array_iterator_get_value(m_tiledb_array_iterator, i,
        reinterpret_cast<const void**>(&field_ptr), &field_size);
    m_cell.set_field_ptr_for_query_idx(i, field_ptr);
    m_cell.set_field_size_in_bytes(i, field_size);
  }
  //Co-ordinates
  tiledb_array_iterator_get_value(m_tiledb_array_iterator, m_num_queried_attributes,
        reinterpret_cast<const void**>(&field_ptr), &field_size);
  assert(field_size == 2u*sizeof(int64_t));
  auto coords_ptr = reinterpret_cast<const int64_t*>(field_ptr);
  m_cell.set_coordinates(coords_ptr[0], coords_ptr[1]);
  return m_cell;
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
  tiledb_array_finalize(tiledb_array);
  if(status == TILEDB_OK)
  {
    auto idx = m_open_arrays_info_vector.size();
    m_open_arrays_info_vector.emplace_back();
    auto& curr_elem = m_open_arrays_info_vector[idx];
    curr_elem.m_idx = idx;
    curr_elem.m_mode = mode_int;
    curr_elem.m_name = array_name;
    curr_elem.m_schema = 0;
    return idx;
  }
  else
    return -1;
}

void VariantStorageManager::close_array(const int ad)
{
  assert(static_cast<size_t>(ad) < m_open_arrays_info_vector.size());
  auto& curr_elem = m_open_arrays_info_vector[ad];
  curr_elem.m_mode = -1;
  curr_elem.m_name.clear();
  curr_elem.m_schema = 0;
}

int VariantStorageManager::define_array(const VariantArraySchema* variant_array_schema)
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
      1000,
      // Compression
      &(compression[0])
  );
  /* Create the array schema */
  return tiledb_array_create(m_tiledb_ctx, &array_schema); 
}

int VariantStorageManager::get_array_schema(const int ad, VariantArraySchema* variant_array_schema)
{
  VERIFY_OR_THROW(static_cast<size_t>(ad) < m_open_arrays_info_vector.size() &&
      m_open_arrays_info_vector[ad].m_name.length());
  auto& curr_elem = m_open_arrays_info_vector[ad];
  // Get array schema
  TileDB_ArraySchema tiledb_array_schema;
  auto status = tiledb_array_load_schema(m_tiledb_ctx, (m_workspace+'/'+curr_elem.m_name).c_str(),
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
  auto dim_domains = std::vector<std::pair<int64_t,int64_t>>(2u*tiledb_array_schema.dim_num_);
  auto dim_domains_ptr = reinterpret_cast<const int64_t*>(tiledb_array_schema.domain_);
  for(auto i=0;i<tiledb_array_schema.dim_num_;++i)
  {
    dim_names[i] = tiledb_array_schema.dimensions_[i];
    dim_domains[i].first = dim_domains_ptr[2*i];
    dim_domains[i].second = dim_domains_ptr[2*i+1];
  }
  *variant_array_schema = std::move(VariantArraySchema(
        curr_elem.m_name,
        attribute_names,
        dim_names,
        dim_domains,
        attribute_types,
        val_num, 
        compression,
        TILEDB_COL_MAJOR)
      );
  curr_elem.m_schema = variant_array_schema;
  // Free array schema
  tiledb_array_free_schema(&tiledb_array_schema);
  return TILEDB_OK;
}

VariantArrayCellIterator VariantStorageManager::begin(
    int ad, const int64_t* range, const std::vector<int>& attribute_ids) const
{
  VERIFY_OR_THROW(static_cast<size_t>(ad) < m_open_arrays_info_vector.size() &&
      m_open_arrays_info_vector[ad].m_name.length());
  auto& curr_elem = m_open_arrays_info_vector[ad];
  assert(curr_elem.m_schema);
  return VariantArrayCellIterator(m_tiledb_ctx, *(curr_elem.m_schema), m_workspace+'/'+curr_elem.m_name,
      range, attribute_ids, m_segment_size);   
}
