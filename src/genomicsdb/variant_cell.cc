#include "variant_cell.h"
#include "variant_field_data.h"
#include "variant_query_config.h"

BufferVariantCell::BufferVariantCell(const VariantArraySchema& array_schema)
{
  clear();
  m_array_schema = &array_schema;
  resize(array_schema.attribute_num());
  for(auto i=0u;i<array_schema.attribute_num();++i)
    update_field_info(i, i);
}

BufferVariantCell::BufferVariantCell(const VariantArraySchema& array_schema, const VariantQueryConfig& query_config)
{
  clear();
  m_array_schema = &array_schema;
  assert(query_config.is_bookkeeping_done());
  resize(query_config.get_num_queried_attributes());
  for(auto i=0u;i<query_config.get_num_queried_attributes();++i)
    update_field_info(i, query_config.get_schema_idx_for_query_idx(i));
}

BufferVariantCell::BufferVariantCell(const VariantArraySchema& array_schema, const std::vector<int>& attribute_ids)
{
  clear();
  m_array_schema = &array_schema;
  resize(attribute_ids.size());
  for(auto i=0u;i<attribute_ids.size();++i)
    update_field_info(i, attribute_ids[i]);
}

void BufferVariantCell::clear()
{
  m_field_ptrs.clear();
  m_field_element_sizes.clear();
  m_field_length_descriptors.clear();
  m_field_lengths.clear();
  m_row_idx = m_begin_column_idx = -1ll;
}

void BufferVariantCell::resize(const size_t num_fields)
{
  //Resize vectors
  m_field_ptrs.resize(num_fields);
  m_field_element_sizes.resize(num_fields);
  m_field_length_descriptors.resize(num_fields);
  m_field_lengths.resize(num_fields);
}
    
void BufferVariantCell::update_field_info(const int query_idx, const int schema_idx)
{
  assert(static_cast<const size_t>(query_idx) < m_field_ptrs.size());
  auto type_index = m_array_schema->type(schema_idx);
  m_field_element_sizes[query_idx] = VariantFieldTypeUtil::size(type_index);
  m_field_length_descriptors[query_idx] = m_array_schema->val_num(schema_idx);
  m_field_lengths[query_idx] = m_field_length_descriptors[query_idx];
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
