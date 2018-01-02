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

#include <algorithm>
#include "variant_query_config.h"
#include "query_variants.h"

using namespace std;

bool ColumnRangeCompare(const ColumnRange& x, const ColumnRange& y)
{
  return (x.first < y.first);
}

void VariantQueryConfig::add_attribute_to_query(const string& name, unsigned schema_idx)
{
  if(m_query_attribute_name_to_query_idx.find(name) == m_query_attribute_name_to_query_idx.end())
  {
    auto idx = m_query_attributes_info_vec.size();
    m_query_attributes_info_vec.emplace_back(name, schema_idx);
    m_query_attribute_name_to_query_idx[name] = idx;
  }
}

void VariantQueryConfig::set_attributes_to_query(const vector<string>& attributeNames)
{
  for(auto i=0u;i<attributeNames.size();++i)
    add_attribute_to_query(attributeNames[i], UNDEFINED_ATTRIBUTE_IDX_VALUE);
}

void VariantQueryConfig::set_rows_to_query(const vector<int64_t>& rowIdxVec)
{
  m_query_rows.resize(rowIdxVec.size());
  for(auto i=0u;i<rowIdxVec.size();++i)
    m_query_rows[i] = rowIdxVec[i];
  std::sort(m_query_rows.begin(), m_query_rows.end());    //useful in querying
  m_query_all_rows = false;
}

void VariantQueryConfig::add_column_interval_to_query(const int64_t colBegin, const int64_t colEnd)
{
  m_query_column_intervals.push_back(make_pair(colBegin, colEnd));
  std::sort(m_query_column_intervals.begin(), m_query_column_intervals.end(), ColumnRangeCompare);
}

void VariantQueryConfig::set_column_interval_to_query(const int64_t colBegin, const int64_t colEnd)
{
  m_query_column_intervals.resize(1u);
  m_query_column_intervals[0] = make_pair(colBegin, colEnd);
}

void VariantQueryConfig::invalidate_array_row_idx_to_query_row_idx_map(bool all_rows)
{
  if(all_rows)
    for(auto i=0ull;i<get_num_rows_in_array();++i)
      m_array_row_idx_to_query_row_idx[i] = UNDEFINED_NUM_ROWS_VALUE;
  else
    for(auto i=0ull;i<get_num_rows_to_query();++i)
    {
      assert(get_array_row_idx_for_query_row_idx(i) >= m_smallest_row_idx
          && (get_array_row_idx_for_query_row_idx(i)-m_smallest_row_idx) < static_cast<int64_t>(m_array_row_idx_to_query_row_idx.size()));
      m_array_row_idx_to_query_row_idx[get_array_row_idx_for_query_row_idx(i)-m_smallest_row_idx] = UNDEFINED_NUM_ROWS_VALUE;
    }
}

void VariantQueryConfig::setup_array_row_idx_to_query_row_idx_map()
{
  if(m_query_all_rows)  //if querying all rows, don't even bother setting up map
    return;
  m_array_row_idx_to_query_row_idx.resize(get_num_rows_in_array());
  invalidate_array_row_idx_to_query_row_idx_map(true);
  //Some queried row idxs may be out of bounds - ignore them
  //This is done by creating a tmp_vector
  std::vector<int64_t> tmp_vector(m_query_rows.size());
  auto valid_rows_idx = 0ull;
  for(auto i=0ull;i<m_query_rows.size();++i)
  {
    //Within bounds
    if(get_array_row_idx_for_query_row_idx(i) >= m_smallest_row_idx &&
        get_array_row_idx_for_query_row_idx(i) < static_cast<int64_t>(get_num_rows_in_array()+m_smallest_row_idx))
    {
      m_array_row_idx_to_query_row_idx[m_query_rows[i]-m_smallest_row_idx] = valid_rows_idx;
      tmp_vector[valid_rows_idx++] = m_query_rows[i];
    }
  }
  tmp_vector.resize(valid_rows_idx);
  m_query_rows = std::move(tmp_vector);
}

void VariantQueryConfig::update_rows_to_query(const std::vector<int64_t>& rows)
{
  assert(is_bookkeeping_done());
  //invalidate old queried rows, if a subset of rows were being queried earlier
  if(!m_query_all_rows)
  {
    invalidate_array_row_idx_to_query_row_idx_map(false); 
    set_rows_to_query(rows);
  }
  else
  {
    set_rows_to_query(rows);
    setup_array_row_idx_to_query_row_idx_map();
  }
  //setup mapping for newly querid rows
  //Some queried row idxs may be out of bounds - ignore them
  //This is done by creating a tmp_vector
  std::vector<int64_t> tmp_vector(m_query_rows.size());
  auto valid_rows_idx = 0ull;
  for(auto i=0ull;i<m_query_rows.size();++i)
  {
    //Within bounds
    if(get_array_row_idx_for_query_row_idx(i) >= m_smallest_row_idx &&
        get_array_row_idx_for_query_row_idx(i) < static_cast<int64_t>(get_num_rows_in_array()+m_smallest_row_idx))
    {
      m_array_row_idx_to_query_row_idx[m_query_rows[i]-m_smallest_row_idx] = valid_rows_idx;
      tmp_vector[valid_rows_idx++] = m_query_rows[i];
    }
  }
  tmp_vector.resize(valid_rows_idx);
  m_query_rows = std::move(tmp_vector);
}

void VariantQueryConfig::update_rows_to_query_to_all_rows()
{
  if(m_query_all_rows)  //already querying all rows
    return;
  assert(is_bookkeeping_done());
  invalidate_array_row_idx_to_query_row_idx_map(false); //invalidate mappings for currently queried rows
  m_query_all_rows = true;
}

void VariantQueryConfig::reorder_query_fields()
{
  auto special_field_names = vector<string>{ "END", "REF", "ALT" };
  m_first_normal_field_query_idx = 0u;
  for(auto i=0u;i<special_field_names.size();++i)
  {
    unsigned query_idx = 0u;
    auto& curr_field_name = special_field_names[i];
    if(get_query_idx_for_name(curr_field_name, query_idx))
    {
      assert(query_idx >= m_first_normal_field_query_idx);
      if(query_idx > m_first_normal_field_query_idx) // == implies already in right place
      {
        auto& other_field_name = m_query_attributes_info_vec[m_first_normal_field_query_idx].m_name;
        //Before swap, update name mappings
        m_query_attribute_name_to_query_idx[curr_field_name] = m_first_normal_field_query_idx;
        m_query_attribute_name_to_query_idx[other_field_name] = query_idx;
        //Swap positions in schema idx and attribute names vector
        std::swap(m_query_attributes_info_vec[query_idx], m_query_attributes_info_vec[m_first_normal_field_query_idx]);
        //Now other_field_name can no longer be used
      }
      ++m_first_normal_field_query_idx;
    }
  }
}

void VariantQueryConfig::flatten_composite_fields(const VidMapper& vid_mapper)
{
  auto has_composite_fields = false;
  for(auto i=0u;i<get_num_queried_attributes();++i)
  {
    auto field_info = vid_mapper.get_field_info(m_query_attributes_info_vec[i].m_name);
    if(field_info == 0)
      throw UnknownQueryAttributeException(std::string("Field ")
          +m_query_attributes_info_vec[i].m_name+" not found in vid mapping");
    auto num_elements_in_tuple = field_info->get_genomicsdb_type().get_num_elements_in_tuple();
    //Multi-element tuple - composite field
    if(num_elements_in_tuple > 1u)
    {
      has_composite_fields = true;
      for(auto j=0u;j<num_elements_in_tuple;++j)
      {
        auto flattened_field_info = vid_mapper.get_flattened_field_info(field_info, j);
        add_attribute_to_query(flattened_field_info->m_name, UNDEFINED_ATTRIBUTE_IDX_VALUE);
      }
    }
  }
  if(has_composite_fields)
  {
    auto num_composite_fields_seen = 0u;
    //Remove composite fields, shift up m_query_attributes_info_vec and re-assign ids
    for(auto i=0u;i<get_num_queried_attributes();++i)
    {
      auto& field_name = m_query_attributes_info_vec[i].m_name;
      //Multi-element tuple - composite field
      if(vid_mapper.get_field_info(field_name)->get_genomicsdb_type().get_num_elements_in_tuple()
          > 1u)
        ++num_composite_fields_seen;
      else
        if(num_composite_fields_seen > 0u) //shift-up elements in m_query_attributes_info_vec
        {
          assert(m_query_attribute_name_to_query_idx[field_name] == i);
          m_query_attribute_name_to_query_idx[field_name] -= num_composite_fields_seen;
          m_query_attributes_info_vec[i-num_composite_fields_seen] =
            std::move(m_query_attributes_info_vec[i]);
        }
    }
    m_query_attributes_info_vec.resize(m_query_attributes_info_vec.size()-num_composite_fields_seen,
        VariantQueryFieldInfo("", UNDEFINED_ATTRIBUTE_IDX_VALUE));
  }
}
