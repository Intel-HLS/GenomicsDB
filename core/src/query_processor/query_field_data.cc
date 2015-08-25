#include "query_field_data.h"
#include "query_config.h"

using namespace std;

//Map from type_index to size_t
unordered_map<type_index, size_t> QueryFieldData::m_type_to_size = unordered_map<type_index, size_t>
{
  { type_index(typeid(int)), sizeof(int) },
  { type_index(typeid(char)), sizeof(char) },
  { type_index(typeid(float)), sizeof(float) },
  { type_index(typeid(double)), sizeof(double) },
  { type_index(typeid(int64_t)), sizeof(int64_t) },
  { type_index(typeid(uint64_t)), sizeof(uint64_t) },
  { type_index(typeid(void)), 0ull }    //opaque type, could be anything
};

QueryFieldData::QueryFieldData()
  : m_element_type(typeid(void))                //undefined
{
  clear();
  m_is_data_ptr_allocated = false;
  reset();
}

void QueryFieldData::reset()
{
  assert(!m_is_data_ptr_allocated);
  m_element_size = UNDEFINED_ATTRIBUTE_IDX_VALUE;
  m_num_elements = 0u;
  m_data_ptr = 0;
}
