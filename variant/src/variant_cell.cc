#include "variant_cell.h"
#include "variant_field_data.h"
#include "variant_query_config.h"

BufferVariantCell::BufferVariantCell(const VariantArraySchema& array_schema, const VariantQueryConfig& query_config)
{
  clear();
  m_array_schema = &array_schema;
  assert(query_config.is_bookkeeping_done());
  //Resize vectors
  m_field_ptrs.resize(query_config.get_num_queried_attributes());
  m_field_element_sizes.resize(query_config.get_num_queried_attributes());
  m_field_length_descriptors.resize(query_config.get_num_queried_attributes());
  m_field_lengths.resize(query_config.get_num_queried_attributes());
  for(auto i=0u;i<query_config.get_num_queried_attributes();++i)
  {
    auto schema_idx = query_config.get_schema_idx_for_query_idx(i);
    auto type_index = m_array_schema->type(schema_idx);
    m_field_element_sizes[i] = VariantFieldTypeUtil::size(type_index);
    m_field_length_descriptors[i] = m_array_schema->val_num(schema_idx);
    m_field_lengths[i] = m_field_length_descriptors[i];
  }
}

BufferVariantCell::BufferVariantCell(const VariantArraySchema& array_schema, const std::vector<int>& attribute_ids)
{
  clear();
  m_array_schema = &array_schema;
  //Resize vectors
  m_field_ptrs.resize(attribute_ids.size());
  m_field_element_sizes.resize(attribute_ids.size());
  m_field_length_descriptors.resize(attribute_ids.size());
  m_field_lengths.resize(attribute_ids.size());
  for(auto i=0u;i<attribute_ids.size();++i)
  {
    auto schema_idx = attribute_ids[i];
    auto type_index = m_array_schema->type(schema_idx);
    m_field_element_sizes[i] = VariantFieldTypeUtil::size(type_index);
    m_field_length_descriptors[i] = m_array_schema->val_num(schema_idx);
    m_field_lengths[i] = m_field_length_descriptors[i];
  }
}

void BufferVariantCell::clear()
{
  m_field_ptrs.clear();
  m_field_element_sizes.clear();
  m_field_length_descriptors.clear();
  m_field_lengths.clear();
  m_row_idx = m_begin_column_idx = -1ll;
}

void BufferVariantCell::set_cell(const void* ptr)
{
  assert(ptr);
  auto cell_ptr = reinterpret_cast<const uint8_t*>(ptr);
  m_row_idx = *(reinterpret_cast<const int64_t*>(cell_ptr));
  m_begin_column_idx = *(reinterpret_cast<const int64_t*>(cell_ptr+sizeof(int64_t))); 
  //Go past co-ordinates and cell size
  uint64_t offset = 2*sizeof(int64_t)+sizeof(size_t);
#ifdef DEBUG
  auto cell_size = *(reinterpret_cast<const size_t*>(cell_ptr+2*sizeof(int64_t)));
#endif
  for(auto i=0u;i<m_field_ptrs.size();++i)
  {
    auto length = m_field_length_descriptors[i];
    //check if variable length field - read length from buffer
    if(length == TILEDB_VAR_NUM)
    {
      length = *(reinterpret_cast<const int*>(cell_ptr+offset));
      m_field_lengths[i] = length;
      offset += sizeof(int);
    }
    m_field_ptrs[i] = cell_ptr + offset;      //field pointer points to region in buffer AFTER the length
    offset += (length*m_field_element_sizes[i]);
  }
#ifdef DEBUG
  assert(offset == cell_size);
#endif
}
