#include <algorithm>
#include "query_config.h"

using namespace std;

bool ColumnRangeCompare(const ColumnRange& x, const ColumnRange& y)
{
  return (x.first < y.first);
}

QueryConfig::QueryConfig() {
  clear();
  m_query_all_rows = true;
}

void QueryConfig::clear() {
  m_query_attributes_names.clear();
  m_query_attributes_schema_idxs.clear();
  m_query_attribute_name_to_query_idx.clear();
  m_query_rows.clear();
  m_query_column_ranges.clear();
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

void QueryConfig::set_rows_to_query(const vector<int64_t>& rowIdxVec)
{
  m_query_rows.resize(rowIdxVec.size());
  for(auto i=0u;i<rowIdxVec.size();++i)
    m_query_rows[i] = rowIdxVec[i];
  std::sort(m_query_rows.begin(), m_query_rows.end());    //useful in querying
  m_query_all_rows = false;
}

void QueryConfig::add_column_range_to_query(const int64_t colBegin, const int64_t colEnd)
{
  m_query_column_ranges.push_back(make_pair(colBegin, colEnd));
  std::sort(m_query_column_ranges.begin(), m_query_column_ranges.end(), ColumnRangeCompare);
}


