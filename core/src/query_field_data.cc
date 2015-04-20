#include "query_field_data.h"
#include "query_config.h"

QueryFieldData::QueryFieldData()
{
  clear();
  m_is_data_ptr_allocated = false;
  reset();
}

void QueryFieldData::reset()
{
  assert(!m_is_data_ptr_allocated);
  m_element_size = UNDEFINED_ATTRIBUTE_IDX_VALUE;
  m_element_type = &(typeid(void));
  m_schema_idx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
  m_num_elements = 0u;
  m_data_ptr = 0;
}
