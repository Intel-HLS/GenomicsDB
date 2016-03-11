#include "variant_cell.h"
#include "variant_field_data.h"
#include "variant_query_config.h"

BufferVariantCell::BufferVariantCell(const VariantArraySchema& array_schema, const VariantQueryConfig& query_config)
{
  clear();
  m_array_schema = &array_schema;
  m_query_config = &query_config;
  m_cell_ptr = 0;
  assert(m_query_config->is_bookkeeping_done());
  //Resize vectors
  m_field_ptrs.resize(m_query_config->get_num_queried_attributes());
  m_field_element_sizes.resize(m_query_config->get_num_queried_attributes());
  m_field_length_descriptors.resize(m_query_config->get_num_queried_attributes());
  m_field_lengths.resize(m_query_config->get_num_queried_attributes());
  for(auto i=0u;i<m_query_config->get_num_queried_attributes();++i)
  {
    auto schema_idx = m_query_config->get_schema_idx_for_query_idx(i);
    auto type_index = std::type_index(*(m_array_schema->type(schema_idx)));
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
}

void BufferVariantCell::set_cell(const void* ptr)
{
  assert(ptr);
  m_cell_ptr = reinterpret_cast<const uint8_t*>(ptr);
  //Go past co-ordinates and cell size
  uint64_t offset = 2*sizeof(int64_t)+sizeof(size_t);
#ifdef DEBUG
  auto cell_size = *(reinterpret_cast<const size_t*>(m_cell_ptr+2*sizeof(int64_t)));
#endif
  for(auto i=0u;i<m_field_ptrs.size();++i)
  {
    auto length = m_field_length_descriptors[i];
    //check if variable length field - read length from buffer
    if(length == VAR_SIZE)
    {
      length = *(reinterpret_cast<const int*>(m_cell_ptr+offset));
      m_field_lengths[i] = length;
      offset += sizeof(int);
    }
    m_field_ptrs[i] = m_cell_ptr + offset;      //field pointer points to region in buffer AFTER the length
    offset += (length*m_field_element_sizes[i]);
  }
#ifdef DEBUG
  assert(offset == cell_size);
#endif
}
