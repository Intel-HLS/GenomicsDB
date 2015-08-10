#include "query_config.h"

using namespace std;

QueryConfig::QueryConfig() {
  clear();
}

void QueryConfig::clear() {
  m_query_attributes_names.clear();
  m_query_attributes_schema_idxs.clear();
  m_query_attribute_name_to_query_idx.clear();
}

void QueryConfig::add_attribute_to_query(const string& name, unsigned schema_idx)
{
  if(m_query_attribute_name_to_query_idx.find(name) == m_query_attribute_name_to_query_idx.end())
  {
    m_query_attributes_names.push_back(name);
    m_query_attributes_schema_idxs.push_back(schema_idx);
    m_query_attribute_name_to_query_idx[name] = m_query_attributes_names.size()-1;
  }
}

void QueryConfig::set_attributes_to_query(const vector<string>& attributeNames)
{
  for(auto i=0u;i<attributeNames.size();++i)
    add_attribute_to_query(attributeNames[i], UNDEFINED_ATTRIBUTE_IDX_VALUE);
}

