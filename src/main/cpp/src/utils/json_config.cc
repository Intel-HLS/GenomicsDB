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

//Enable asserts
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <zlib.h>
#include "json_config.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw RunConfigException(#X);

const char* g_json_indent_unit = "    ";

void JSONConfigBase::clear()
{
  m_workspaces.clear();
  m_array_names.clear();
  m_column_ranges.clear();
  m_row_ranges.clear();
  m_attributes.clear();
  m_sorted_column_partitions.clear();
  m_sorted_row_partitions.clear();
  m_vid_mapping_file.clear();
  m_callset_mapping_file.clear();
}

ColumnRange verify_contig_position_and_get_tiledb_column_interval(const ContigInfo& contig_info,
    const int64_t begin, const int64_t end)
{
  ColumnRange result;
  result.first = std::min(begin, end);
  result.second = std::max(begin, end);
  if(result.first > contig_info.m_length)
    throw RunConfigException(std::string("Position ")+std::to_string(result.first)
        +" queried for contig "+contig_info.m_name+" which is of length "
        +std::to_string(contig_info.m_length)+"; queried position is past end of contig");
  if(result.second > contig_info.m_length)
  {
    std::cerr << "WARNING: position "+std::to_string(result.second)
        +" queried for contig "+contig_info.m_name+" which is of length "
        +std::to_string(contig_info.m_length)+"; queried interval is past end of contig, truncating to contig length";
    result.second = contig_info.m_length;
  }
  // Subtract 1 as TileDB is 0-based and genomics (VCF) is 1-based
  result.first += (contig_info.m_tiledb_column_offset - 1);
  result.second += (contig_info.m_tiledb_column_offset - 1);
  return result;
}

void JSONConfigBase::extract_contig_interval_from_object(const rapidjson::Value& curr_json_object,
    const VidMapper* id_mapper, ColumnRange& result)
{
  // This MUST be a dictionary of the form "contig" : position or "contig" : [start, end]
  VERIFY_OR_THROW(curr_json_object.IsObject());
  VERIFY_OR_THROW(curr_json_object.MemberCount() == 1);
  auto itr = curr_json_object.MemberBegin();
  std::string contig_name = itr->name.GetString();
  const rapidjson::Value* contig_position_ptr = 0;
  if(contig_name == "contig_position") // Produced by PB
  {
    VERIFY_OR_THROW(itr->value.IsObject());
    auto& pb_contig_position_dict = itr->value;
    VERIFY_OR_THROW(pb_contig_position_dict.MemberCount() == 2);
    VERIFY_OR_THROW(pb_contig_position_dict.HasMember("contig")
        && pb_contig_position_dict["contig"].IsString()
        && pb_contig_position_dict.HasMember("position")
        && pb_contig_position_dict["position"].IsInt64());
    contig_name = pb_contig_position_dict["contig"].GetString();
    contig_position_ptr = &(pb_contig_position_dict["position"]);
  }
  else
    contig_position_ptr = &(itr->value);
  ContigInfo contig_info;
  VERIFY_OR_THROW(id_mapper != 0 && id_mapper->is_initialized());
  if (!id_mapper->get_contig_info(contig_name, contig_info))
    throw VidMapperException("JSONConfigBase::read_from_file: Invalid contig name : " + contig_name );
  const rapidjson::Value& contig_position = *contig_position_ptr;
  if (contig_position.IsArray())
  {
    VERIFY_OR_THROW(contig_position.Size() == 2);
    VERIFY_OR_THROW(contig_position[0u].IsInt64());
    VERIFY_OR_THROW(contig_position[1u].IsInt64());
    result = verify_contig_position_and_get_tiledb_column_interval(contig_info,
        contig_position[0u].GetInt64(), contig_position[1u].GetInt64());
  }
  else // single position
  {
    VERIFY_OR_THROW(contig_position.IsInt64());
    auto contig_position_int = contig_position.GetInt64();
    result = verify_contig_position_and_get_tiledb_column_interval(contig_info,
        contig_position_int, contig_position_int);
  }
}

//JSON produced by Protobuf - { "low":<>, "high":<> }
bool JSONConfigBase::extract_interval_from_PB_struct_or_return_false(const rapidjson::Value& curr_json_object,
    const VidMapper* id_mapper,
    ColumnRange& result)
{
  VERIFY_OR_THROW(id_mapper != 0 && id_mapper->is_initialized());
  if(curr_json_object.IsObject())
  {
    //Dictionary of the form { "high": <>, "low": <> }
    if(curr_json_object.MemberCount() == 2u &&
        curr_json_object.HasMember("low") && curr_json_object.HasMember("high"))
    {
      result.first = curr_json_object["low"].GetInt64();
      result.second = curr_json_object["high"].GetInt64();
      return true;
    }
    if(curr_json_object.MemberCount() == 1u)
    {
      //Dictionary representing column interval
      if(curr_json_object.HasMember("column_interval"))
      {
        const auto& interval_object = curr_json_object["column_interval"];
        //This could be TileDB column interval or contig interval
        if(interval_object.IsObject())
        {
          //TileDB column interval
          if(interval_object.HasMember("column_interval")
              && interval_object["column_interval"].IsObject()
              && interval_object["column_interval"].HasMember("begin")
              && interval_object["column_interval"].HasMember("end"))
          {
            result.first = interval_object["column_interval"]["begin"].GetInt64();
            result.second = interval_object["column_interval"]["end"].GetInt64();
            return true;
          }
          else
            if(interval_object.HasMember("contig_interval")
                && interval_object["contig_interval"].IsObject()
                && interval_object["contig_interval"].HasMember("contig")
                && interval_object["contig_interval"].HasMember("begin")
                && interval_object["contig_interval"].HasMember("end"))
            {
              ContigInfo contig_info;
              auto contig_name = interval_object["contig_interval"]["contig"].GetString();
              if (!id_mapper->get_contig_info(contig_name, contig_info))
                throw VidMapperException(std::string("JSONConfigBase::read_from_file: Invalid contig name : ")
                    + contig_name );
              result = verify_contig_position_and_get_tiledb_column_interval(contig_info,
                  interval_object["contig_interval"]["begin"].GetInt64(),
                  interval_object["contig_interval"]["end"].GetInt64());
              return true;
            }
        }
      }
      else
      {
        //Dictionary representing column position
        if(curr_json_object.HasMember("column"))
        {
          const auto& interval_object = curr_json_object["column"];
          //This could be TileDB column or contig position
          if(interval_object.IsObject())
          {
            //TileDB column
            if(interval_object.HasMember("tiledb_column")
                && interval_object["tiledb_column"].IsInt64())
            {
              result.first = interval_object["tiledb_column"].GetInt64();
              result.second = result.first;
              return true;
            }
            else
              if(interval_object.HasMember("contig_position")
                  && interval_object["contig_position"].IsObject()
                  && interval_object["contig_position"].HasMember("contig")
                  && interval_object["contig_position"].HasMember("position"))
              {
                ContigInfo contig_info;
                auto contig_name = interval_object["contig_position"]["contig"].GetString();
                if (!id_mapper->get_contig_info(contig_name, contig_info))
                  throw VidMapperException(std::string("JSONConfigBase::read_from_file: Invalid contig name : ")
                      + contig_name );
                auto contig_position_int = interval_object["contig_position"]["position"].GetInt64();
                result = verify_contig_position_and_get_tiledb_column_interval(contig_info,
                    contig_position_int, contig_position_int);
                return true;
              }
          }
        }
      }
    }
  }
  return false;
}

void JSONConfigBase::read_from_file(const std::string& filename, const VidMapper* id_mapper, const int rank)
{
  std::ifstream ifs(filename.c_str());
  VERIFY_OR_THROW(ifs.is_open());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  m_json.Parse(str.c_str());
  if(m_json.HasParseError())
    throw RunConfigException(std::string("Syntax error in JSON file ")+filename);
  VidMapper tmp_vid_mapper;
  //Null or un-initialized
  if(id_mapper == 0 || !(id_mapper->is_initialized()))
  {
    read_and_initialize_vid_and_callset_mapping_if_available(&tmp_vid_mapper, rank);
    if(tmp_vid_mapper.is_initialized())
      id_mapper = &tmp_vid_mapper;
  }
  //Workspace
  if(m_json.HasMember("workspace"))
  {
    const rapidjson::Value& workspace = m_json["workspace"];
    //workspace could be an array, one workspace dir for every rank
    if(workspace.IsArray())
    {
      for(rapidjson::SizeType i=0;i<workspace.Size();++i)
      {
        VERIFY_OR_THROW(workspace[i].IsString());
        m_workspaces.push_back(workspace[i].GetString());
      }
    }
    else //workspace is simply a string
    {
      VERIFY_OR_THROW(workspace.IsString());
      m_workspaces.push_back(workspace.GetString());
      m_single_workspace_path = true;
    }
  }
  //Array
  VERIFY_OR_THROW(!(m_json.HasMember("array") && m_json.HasMember("array_name")));
  if(m_json.HasMember("array") || m_json.HasMember("array_name"))
  {
    const rapidjson::Value& array_name = m_json.HasMember("array") ? m_json["array"] : m_json["array_name"];
    //array could be an array, one array dir for every rank
    if(array_name.IsArray())
    {
      for(rapidjson::SizeType i=0;i<array_name.Size();++i)
      {
        VERIFY_OR_THROW(array_name[i].IsString());
        m_array_names.push_back(array_name[i].GetString());
      }
    }
    else //array is simply a string
    {
      VERIFY_OR_THROW(array_name.IsString());
      m_array_names.push_back(array_name.GetString());
      m_single_array_name = true;
    }
  }
  //VERIFY_OR_THROW(m_json.HasMember("query_column_ranges") || m_json.HasMember("column_partitions")
  //|| m_json.HasMember("scan_full"));
  if(m_json.HasMember("scan_full"))
  {
    m_scan_whole_array = true;
    m_column_ranges.resize(1u, std::vector<ColumnRange>({ColumnRange(0, INT64_MAX-1)}));
    m_single_query_column_ranges_vector = true;
    m_row_ranges.resize(1u, std::vector<RowRange>({RowRange(0, INT64_MAX-1)}));
    m_single_query_row_ranges_vector = true;
  }
  else
  {
    VERIFY_OR_THROW(!(m_json.HasMember("row_partitions") && m_json.HasMember("column_partitions"))
        && "Cannot have both \"row_partitions\" and \"column_partitions\" simultaneously in the JSON file");
    VERIFY_OR_THROW((m_json.HasMember("query_column_ranges") || m_json.HasMember("column_partitions") ||
          m_json.HasMember("query_row_ranges") || m_json.HasMember("row_partitions")) &&
        "Must have one of \"query_column_ranges\" or \"column_partitions\" or \"query_row_ranges\" or \"row_partitions\"");
    VERIFY_OR_THROW((!m_json.HasMember("query_column_ranges") || !m_json.HasMember("column_partitions")) &&
        "Cannot use both \"query_column_ranges\" and \"column_partitions\" simultaneously");
    //Query columns
    //Example:  [ [ [0,5], 45 ], [ 76, 87 ] ]
    //This means that rank 0 will have 2 query intervals: [0-5] and [45-45] and rank 1 will have
    //2 intervals [76-76] and [87-87]
    //But you could have a single innermost list - with this option all ranks will query the same list 
    if(m_json.HasMember("query_column_ranges"))
    {
      const rapidjson::Value& q1 = m_json["query_column_ranges"];
      VERIFY_OR_THROW(q1.IsArray());
      if(q1.Size() == 1)
        m_single_query_column_ranges_vector = true;
      m_column_ranges.resize(q1.Size());
      for(rapidjson::SizeType i=0;i<q1.Size();++i)
      {
        auto& curr_q1_entry = q1[i];
        VERIFY_OR_THROW(curr_q1_entry.IsArray() || curr_q1_entry.IsObject());
        //JSON produced by Protobuf - { "range_list": [ { "low":<>, "high":<> } ] }
        if(curr_q1_entry.IsObject())
        {
          if(curr_q1_entry.MemberCount() == 0u)
            continue;
          VERIFY_OR_THROW(curr_q1_entry.MemberCount() == 1u &&
              (curr_q1_entry.HasMember("range_list")
               || curr_q1_entry.HasMember("column_or_interval_list")));
        }
        const rapidjson::Value& q2 = curr_q1_entry.IsArray() ? q1[i]
          : (curr_q1_entry.HasMember("range_list")
              ?  curr_q1_entry["range_list"]
              : curr_q1_entry["column_or_interval_list"]);
        VERIFY_OR_THROW(q2.IsArray());
        m_column_ranges[i].resize(q2.Size());
        for(rapidjson::SizeType j=0;j<q2.Size();++j)
        {
          const rapidjson::Value& q3 = q2[j];
          //q3 is list of 2 elements to represent query interval
          if(q3.IsArray())
          {
            VERIFY_OR_THROW(q3.Size() == 2);
            VERIFY_OR_THROW(q3[0u].IsInt64());
            VERIFY_OR_THROW(q3[1u].IsInt64());
            m_column_ranges[i][j].first = q3[0u].GetInt64();
            m_column_ranges[i][j].second = q3[1u].GetInt64();
          }
          else if (q3.IsInt64()) //single position in Tile DB format
          {
            m_column_ranges[i][j].first = q3.GetInt64();
            m_column_ranges[i][j].second = q3.GetInt64();
          }
          else if (q3.IsString()) // Query is for entire contig
          {
            ContigInfo contig_info;
            std::string contig_name = q3.GetString();
            assert(id_mapper != 0);
            if (!id_mapper->get_contig_info(contig_name, contig_info))
              throw VidMapperException("JSONConfigBase::read_from_file: Invalid contig name : " + contig_name );
            m_column_ranges[i][j].first = contig_info.m_tiledb_column_offset;
            m_column_ranges[i][j].second = contig_info.m_tiledb_column_offset + contig_info.m_length - 1;
          }
          else
            if(!extract_interval_from_PB_struct_or_return_false(q3, id_mapper, m_column_ranges[i][j])) //check if PB based JSON
            {
              //must be object { "chr" : [ b , e ] }
              extract_contig_interval_from_object(q3, id_mapper, m_column_ranges[i][j]);
            }
          if(m_column_ranges[i][j].first > m_column_ranges[i][j].second)
            std::swap<int64_t>(m_column_ranges[i][j].first, m_column_ranges[i][j].second);
        }
      }
    }
    else
      if (m_json.HasMember("column_partitions"))
      {
        m_column_partitions_specified = true;
        //column_partitions_array itself is an array of the form [ { "begin" : <value> }, { "begin":<value>} ]
        auto& column_partitions_array = m_json["column_partitions"];
        VERIFY_OR_THROW(column_partitions_array.IsArray());
        m_sorted_column_partitions.resize(column_partitions_array.Size());
        m_column_ranges.resize(column_partitions_array.Size());
        std::unordered_map<int64_t, unsigned> begin_to_idx;
        auto workspace_string = m_single_workspace_path ? m_workspaces[0] : "";
        auto array_name_string = m_single_array_name ? m_array_names[0] : "";
        for(rapidjson::SizeType partition_idx=0;partition_idx<column_partitions_array.Size();++partition_idx)
        {
          m_column_ranges[partition_idx].resize(1);      //only 1 std::pair
          //{ "begin" : <Val> }
          const auto& curr_partition_info_dict = column_partitions_array[partition_idx];
          VERIFY_OR_THROW(curr_partition_info_dict.IsObject());
          VERIFY_OR_THROW(curr_partition_info_dict.HasMember(JSON_LOADER_PARTITION_INFO_BEGIN_FIELD_NAME));
          auto& begin_json_value = curr_partition_info_dict[JSON_LOADER_PARTITION_INFO_BEGIN_FIELD_NAME];
          //Either a TileDB column idx or a dictionary of the form { "chr1" : [ 5, 6 ] }
          VERIFY_OR_THROW(begin_json_value.IsInt64() || begin_json_value.IsObject());
          if(begin_json_value.IsInt64())
            m_column_ranges[partition_idx][0].first = begin_json_value.GetInt64();
          else
            extract_contig_interval_from_object(begin_json_value, id_mapper, m_column_ranges[partition_idx][0]);
          m_column_ranges[partition_idx][0].second = INT64_MAX-1;
          if(curr_partition_info_dict.HasMember(JSON_LOADER_PARTITION_INFO_END_FIELD_NAME))
          {
            auto& end_json_value = curr_partition_info_dict[JSON_LOADER_PARTITION_INFO_END_FIELD_NAME];
            //Either a TileDB column idx or a dictionary of the form { "chr1" : [ 5, 6 ] }
            VERIFY_OR_THROW(end_json_value.IsInt64() || end_json_value.IsObject());
            if(end_json_value.IsInt64())
              m_column_ranges[partition_idx][0].second = end_json_value.GetInt64();
            else
            {
              ColumnRange tmp_range;
              extract_contig_interval_from_object(end_json_value, id_mapper, tmp_range);
              m_column_ranges[partition_idx][0].second = tmp_range.first;
            }
          }
          if(m_column_ranges[partition_idx][0].first > m_column_ranges[partition_idx][0].second)
            std::swap<int64_t>(m_column_ranges[partition_idx][0].first, m_column_ranges[partition_idx][0].second);
          if(curr_partition_info_dict.HasMember("workspace"))
          {
            if(column_partitions_array.Size() > m_workspaces.size())
              m_workspaces.resize(column_partitions_array.Size(), workspace_string);
            m_workspaces[partition_idx] = curr_partition_info_dict["workspace"].GetString();
            m_single_workspace_path = false;
          }
          if(curr_partition_info_dict.HasMember("array") || curr_partition_info_dict.HasMember("array_name"))
          {
            VERIFY_OR_THROW(!(curr_partition_info_dict.HasMember("array")
                  && curr_partition_info_dict.HasMember("array_name")));
            if(column_partitions_array.Size() >= m_array_names.size())
              m_array_names.resize(column_partitions_array.Size(), array_name_string);
            m_array_names[partition_idx] = curr_partition_info_dict.HasMember("array")
              ? curr_partition_info_dict["array"].GetString()
              : curr_partition_info_dict["array_name"].GetString();
            m_single_array_name = false;
          }
          //Mapping from begin pos to index
          begin_to_idx[m_column_ranges[partition_idx][0].first] = partition_idx;
          m_sorted_column_partitions[partition_idx].first = m_column_ranges[partition_idx][0].first;
          m_sorted_column_partitions[partition_idx].second = m_column_ranges[partition_idx][0].second;
        }
        //Sort in ascending order
        std::sort(m_sorted_column_partitions.begin(), m_sorted_column_partitions.end(), ColumnRangeCompare);
        //Set end value if not valid
        for(auto i=0ull;i+1u<m_sorted_column_partitions.size();++i)
        {
          VERIFY_OR_THROW(m_sorted_column_partitions[i].first != m_sorted_column_partitions[i+1u].first
              && "Cannot have two column partitions with the same begin value");
          if(m_sorted_column_partitions[i].second >= m_sorted_column_partitions[i+1u].first)
            m_sorted_column_partitions[i].second = m_sorted_column_partitions[i+1u].first-1;
          auto idx = begin_to_idx[m_sorted_column_partitions[i].first];
          m_column_ranges[idx][0].second = m_sorted_column_partitions[i].second;
        }
      }
    VERIFY_OR_THROW((!m_json.HasMember("query_row_ranges") || !m_json.HasMember("row_partitions")) &&
        "Cannot use both \"query_row_ranges\" and \"row_partitions\" simultaneously");
    //Query rows
    //Example:  [ [ [0,5], 45 ], [ 76, 87 ] ]
    //This means that rank 0 will query rows: [0-5] and [45-45] and rank 1 will have
    //2 intervals [76-76] and [87-87]
    //But you could have a single innermost list - with this option all ranks will query the same list 
    if(m_json.HasMember("query_row_ranges"))
    {
      const rapidjson::Value& q1 = m_json["query_row_ranges"];
      VERIFY_OR_THROW(q1.IsArray());
      if(q1.Size() == 1)
        m_single_query_row_ranges_vector = true;
      m_row_ranges.resize(q1.Size());
      for(rapidjson::SizeType i=0;i<q1.Size();++i)
      {
        auto& curr_q1_entry = q1[i];
        VERIFY_OR_THROW(curr_q1_entry.IsArray() || curr_q1_entry.IsObject());
        //JSON produced by Protobuf - { "range_list": [ { "low":<>, "high":<> } ] }
        if(curr_q1_entry.IsObject())
          VERIFY_OR_THROW(curr_q1_entry.MemberCount() == 1u && curr_q1_entry.HasMember("range_list"));
        const rapidjson::Value& q2 = curr_q1_entry.IsArray() ? q1[i]
          : curr_q1_entry["range_list"];
        VERIFY_OR_THROW(q2.IsArray());
        m_row_ranges[i].resize(q2.Size());
        for(rapidjson::SizeType j=0;j<q2.Size();++j)
        {
          const rapidjson::Value& q3 = q2[j];
          //q3 is list of 2 elements to represent query row interval
          if(q3.IsArray())
          {
            VERIFY_OR_THROW(q3.Size() == 2);
            VERIFY_OR_THROW(q3[0u].IsInt64());
            VERIFY_OR_THROW(q3[1u].IsInt64());
            m_row_ranges[i][j].first = q3[0u].GetInt64();
            m_row_ranges[i][j].second = q3[1u].GetInt64();
          }
          else
          {
            if(q3.IsInt64()) //single position
            {
              m_row_ranges[i][j].first = q3.GetInt64();
              m_row_ranges[i][j].second = q3.GetInt64();
            }
            else
              if(q3.IsObject()) //Must be PB generated JSON object { "low": <>, "high": <> }
                VERIFY_OR_THROW(extract_interval_from_PB_struct_or_return_false(q3, id_mapper, m_row_ranges[i][j]));
          }
          if(m_row_ranges[i][j].first > m_row_ranges[i][j].second)
            std::swap<int64_t>(m_row_ranges[i][j].first, m_row_ranges[i][j].second);
        }
      }
    }
    else
      if(m_json.HasMember("row_partitions"))
      {
        m_row_partitions_specified = true;
        //row_partitions value itself is an array [ { "begin" : <value> } ]
        auto& row_partitions_array = m_json["row_partitions"];
        VERIFY_OR_THROW(row_partitions_array.IsArray());
        m_sorted_row_partitions.resize(row_partitions_array.Size());
        m_row_ranges.resize(row_partitions_array.Size());
        std::unordered_map<int64_t, unsigned> begin_to_idx;
        auto workspace_string = m_single_workspace_path ? m_workspaces[0] : "";
        auto array_name_string = m_single_array_name ? m_array_names[0] : "";
        for(rapidjson::SizeType partition_idx=0u;partition_idx<row_partitions_array.Size();++partition_idx)
        {
          const auto& curr_partition_info_dict = row_partitions_array[partition_idx];
          VERIFY_OR_THROW(curr_partition_info_dict.IsObject());
          VERIFY_OR_THROW(curr_partition_info_dict.HasMember(JSON_LOADER_PARTITION_INFO_BEGIN_FIELD_NAME));
          m_row_ranges[partition_idx].resize(1);      //only 1 std::pair
          m_row_ranges[partition_idx][0].first = curr_partition_info_dict[JSON_LOADER_PARTITION_INFO_BEGIN_FIELD_NAME].GetInt64();
          m_row_ranges[partition_idx][0].second = INT64_MAX-1;
          if(curr_partition_info_dict.HasMember(JSON_LOADER_PARTITION_INFO_END_FIELD_NAME))
            m_row_ranges[partition_idx][0].second = curr_partition_info_dict[JSON_LOADER_PARTITION_INFO_END_FIELD_NAME].GetInt64();
          if(m_row_ranges[partition_idx][0].first > m_row_ranges[partition_idx][0].second)
            std::swap<int64_t>(m_row_ranges[partition_idx][0].first, m_row_ranges[partition_idx][0].second);
          if(curr_partition_info_dict.HasMember("workspace"))
          {
            if(row_partitions_array.Size() > m_workspaces.size())
              m_workspaces.resize(row_partitions_array.Size(), workspace_string);
            m_workspaces[partition_idx] = curr_partition_info_dict["workspace"].GetString();
            m_single_workspace_path = false;
          }
          if(curr_partition_info_dict.HasMember("array") || curr_partition_info_dict.HasMember("array_name"))
          {
            VERIFY_OR_THROW(!(curr_partition_info_dict.HasMember("array")
                  && curr_partition_info_dict.HasMember("array_name")));
            if(row_partitions_array.Size() >= m_array_names.size())
              m_array_names.resize(row_partitions_array.Size(), array_name_string);
            m_array_names[partition_idx] = curr_partition_info_dict.HasMember("array")
              ? curr_partition_info_dict["array"].GetString()
              : curr_partition_info_dict["array_name"].GetString();
            m_single_array_name = false;
          }
          //Mapping from begin pos to index
          begin_to_idx[m_row_ranges[partition_idx][0].first] = partition_idx;
          m_sorted_row_partitions[partition_idx].first = m_row_ranges[partition_idx][0].first;
          m_sorted_row_partitions[partition_idx].second = m_row_ranges[partition_idx][0].second;
        }
        //Sort in ascending order
        std::sort(m_sorted_row_partitions.begin(), m_sorted_row_partitions.end(), ColumnRangeCompare);
        //Set end value if not valid
        for(auto i=0ull;i+1u<m_sorted_row_partitions.size();++i)
        {
          VERIFY_OR_THROW(m_sorted_row_partitions[i].first != m_sorted_row_partitions[i+1u].first
              && "Cannot have two row partitions with the same begin value");
          if(m_sorted_row_partitions[i].second >= m_sorted_row_partitions[i+1u].first)
            m_sorted_row_partitions[i].second = m_sorted_row_partitions[i+1u].first-1;
          auto idx = begin_to_idx[m_sorted_row_partitions[i].first];
          m_row_ranges[idx][0].second = m_sorted_row_partitions[i].second;
        }
      }
  }
  if(m_json.HasMember("query_attributes") && m_json.HasMember("attributes"))
    throw RunConfigException("Query configuration cannot have both \"query_attributes\" and \"attributes\"");
  if(m_json.HasMember("query_attributes") || m_json.HasMember("attributes"))
  {
    const rapidjson::Value& q1 = m_json.HasMember("query_attributes") ?
      m_json["query_attributes"] : m_json["attributes"];
    VERIFY_OR_THROW(q1.IsArray());
    m_attributes.resize(q1.Size());
    for(rapidjson::SizeType i=0;i<q1.Size();++i)
    {
      const rapidjson::Value& q2 = q1[i];
      VERIFY_OR_THROW(q2.IsString());
      m_attributes[i] = std::move(std::string(q2.GetString()));
    }
  }
  if(m_json.HasMember("segment_size"))
  {
    VERIFY_OR_THROW(m_json["segment_size"].IsInt64());
    m_segment_size = m_json["segment_size"].GetInt64();
  }
}

const std::string& JSONConfigBase::get_workspace(const int rank) const
{
  VERIFY_OR_THROW((m_single_workspace_path || static_cast<size_t>(rank) < m_workspaces.size())
     && ("Workspace not defined for rank "+std::to_string(rank)).c_str());
  if(m_single_workspace_path)
    return m_workspaces[0];
  return m_workspaces[rank];
}

const std::string& JSONConfigBase::get_array_name(const int rank) const
{
  VERIFY_OR_THROW((m_single_array_name || static_cast<size_t>(rank) < m_array_names.size())
      && ("Could not find array for rank "+std::to_string(rank)).c_str());
  if(m_single_array_name)
    return m_array_names[0];
  return m_array_names[rank];
}

RowRange JSONConfigBase::get_row_partition(const int rank, const unsigned idx) const
{
  if(!m_row_partitions_specified)
    return RowRange(0, INT64_MAX-1);
  auto fixed_rank = m_single_query_row_ranges_vector ? 0 : rank;
  if(static_cast<size_t>(fixed_rank) >= m_row_ranges.size())
    throw RunConfigException(std::string("No row partition/query interval available for process with rank ")
        +std::to_string(rank));
  VERIFY_OR_THROW(idx < m_row_ranges[fixed_rank].size());
  return m_row_ranges[fixed_rank][idx];
}

ColumnRange JSONConfigBase::get_column_partition(const int rank, const unsigned idx) const
{
  if(!m_column_partitions_specified)
    return ColumnRange(0, INT64_MAX-1);
  auto fixed_rank = m_single_query_column_ranges_vector ? 0 : rank;
  if(static_cast<size_t>(fixed_rank) >= m_column_ranges.size())
    throw RunConfigException(std::string("No column partition/query interval available for process with rank ")
        +std::to_string(rank));
  VERIFY_OR_THROW(idx < m_column_ranges[fixed_rank].size());
  return m_column_ranges[fixed_rank][idx];
}

const std::vector<RowRange>& JSONConfigBase::get_query_row_ranges(const int rank) const
{
  auto fixed_rank = m_single_query_row_ranges_vector ? 0 : rank;
  if(static_cast<size_t>(fixed_rank) >= m_row_ranges.size())
    throw RunConfigException(std::string("No row partition/query row range available for process with rank ")
        +std::to_string(rank));
  return m_row_ranges[fixed_rank];
}

const std::vector<ColumnRange>& JSONConfigBase::get_query_column_ranges(const int rank) const
{
  auto fixed_rank = m_single_query_column_ranges_vector ? 0 : rank;
  if(static_cast<size_t>(fixed_rank) >= m_column_ranges.size())
    throw RunConfigException(std::string("No column partition/query column range available for process with rank ")
        +std::to_string(rank));
  return m_column_ranges[fixed_rank];
}

void JSONConfigBase::read_and_initialize_vid_and_callset_mapping_if_available(VidMapper* id_mapper, const int rank)
{
  auto& vid_mapping_file = m_vid_mapping_file;
  auto& callset_mapping_file = m_callset_mapping_file;
  //Callset mapping file and vid file
  if(m_json.HasMember("vid_mapping_file"))
  {
    //Over-ride callset mapping file in top-level config if necessary
    if(m_json.HasMember("callset_mapping_file"))
    {
      const rapidjson::Value& v = m_json["callset_mapping_file"];
      //Could be array - one for each process
      if(v.IsArray())
      {
        VERIFY_OR_THROW(rank < static_cast<int>(v.Size()));
        VERIFY_OR_THROW(v[rank].IsString());
        callset_mapping_file = v[rank].GetString();
      }
      else
      {
        VERIFY_OR_THROW(v.IsString());
        callset_mapping_file = v.GetString();
      }
    }
    const rapidjson::Value& v = m_json["vid_mapping_file"];
    //Could be array - one for each process
    if(v.IsArray())
    {
      VERIFY_OR_THROW(rank < static_cast<int>(v.Size()));
      VERIFY_OR_THROW(v[rank].IsString());
      vid_mapping_file = v[rank].GetString();
    }
    else //or single string for all processes
    {
      VERIFY_OR_THROW(v.IsString());
      vid_mapping_file = v.GetString();
    }
  }
  if(id_mapper && !(vid_mapping_file.empty()))
    (*id_mapper) = std::move(FileBasedVidMapper(vid_mapping_file, callset_mapping_file, m_lb_callset_row_idx, m_ub_callset_row_idx, false));
}

void JSONBasicQueryConfig::update_from_loader(JSONLoaderConfig* loader_config, const int rank)
{
  if (!loader_config)
    return;
  // Update workspace if it is not provided in query json
  if (!m_json.HasMember("workspace"))
  {
    // Set as single workspace
    m_single_workspace_path = true;
    m_workspaces.push_back(loader_config->get_workspace(rank));
  }
  // Update array if it is not provided in query json
  if (!(m_json.HasMember("array") || m_json.HasMember("array_name")))
  {
    // Set as single array
    m_single_array_name = true;
    m_array_names.push_back(loader_config->get_array_name(rank));
  }
  //Vid mapping
  if(!m_json.HasMember("vid_mapping_file"))
    m_vid_mapping_file = loader_config->get_vid_mapping_filename();
  //Callset mapping
  if(!m_json.HasMember("callset_mapping_file"))
    m_callset_mapping_file = loader_config->get_callset_mapping_filename();
}

void JSONBasicQueryConfig::subset_query_column_ranges_based_on_partition(const JSONLoaderConfig* loader_config, const int rank)
{
  if(loader_config)
  {
    // Check if the partitioning is column-based, if so, pick the column corresponding to the rank
    // and update the m_column_ranges
    if (loader_config->is_partitioned_by_column())
    {
      ColumnRange my_rank_loader_column_range = loader_config->get_column_partition(rank);
      std::vector<ColumnRange> my_rank_queried_columns;
      for(auto queried_column_range : get_query_column_ranges(rank))
      {
        if(queried_column_range.second >= my_rank_loader_column_range.first
            && queried_column_range.first <= my_rank_loader_column_range.second)
          my_rank_queried_columns.emplace_back(queried_column_range);
      }
      auto idx = m_single_query_column_ranges_vector ? 0 : rank;
      assert(static_cast<size_t>(idx) < m_column_ranges.size());
      m_column_ranges[idx] = std::move(my_rank_queried_columns);
    }
  }
}

void JSONBasicQueryConfig::read_from_file(const std::string& filename, VariantQueryConfig& query_config,
    VidMapper* id_mapper, const int rank, JSONLoaderConfig* loader_config)
{
  //Need to parse here first because id_mapper initialization in read_and_initialize_vid_and_callset_mapping_if_available() requires
  //valid m_json object
  std::ifstream ifs(filename.c_str());
  VERIFY_OR_THROW(ifs.is_open());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  m_json.Parse(str.c_str());
  if(m_json.HasParseError())
    throw RunConfigException(std::string("Syntax error in JSON file ")+filename);
  // Update from loader_config_file
  update_from_loader(loader_config, rank);
  if (id_mapper)
  {
    read_and_initialize_vid_and_callset_mapping_if_available(id_mapper, rank);
    VERIFY_OR_THROW(id_mapper->is_initialized() && "No valid vid_mapping_file provided in loader or query JSON");
  }
  JSONConfigBase::read_from_file(filename, id_mapper, rank);
  //Workspace
  VERIFY_OR_THROW(m_workspaces.size() && "No workspace specified");
  VERIFY_OR_THROW((m_single_workspace_path || static_cast<size_t>(rank) < m_workspaces.size())
      && ("Could not find workspace for rank "+std::to_string(rank)).c_str());
  auto& workspace = m_single_workspace_path ? m_workspaces[0] : m_workspaces[rank];
  VERIFY_OR_THROW(workspace != "" && "Empty workspace string");
  //Array
  VERIFY_OR_THROW(m_array_names.size() && "No array specified");
  VERIFY_OR_THROW((m_single_array_name || static_cast<size_t>(rank) < m_array_names.size())
      && ("Could not find array for rank "+std::to_string(rank)).c_str());
  auto& array_name = m_single_array_name ? m_array_names[0] : m_array_names[rank];
  VERIFY_OR_THROW(array_name != "" && "Empty array name");
  //Query columns
  VERIFY_OR_THROW((m_column_ranges.size() || m_scan_whole_array) && "Query column ranges not specified");
  subset_query_column_ranges_based_on_partition(loader_config, rank);
  if(!m_scan_whole_array)
  {
    VERIFY_OR_THROW((m_single_query_column_ranges_vector || static_cast<size_t>(rank) < m_column_ranges.size())
        && "Rank >= query column ranges vector size");
    for(const auto& range : get_query_column_ranges(rank))
      query_config.add_column_interval_to_query(range.first, range.second);
  }
  //Query rows
  if(!m_scan_whole_array && m_row_ranges.size())
  {
    VERIFY_OR_THROW((m_single_query_row_ranges_vector || static_cast<size_t>(rank) < m_row_ranges.size())
        && "Rank >= query row ranges vector size");
    std::vector<int64_t> row_idxs;
    for(const auto& range : get_query_row_ranges(rank))
    {
      auto j = row_idxs.size();
      row_idxs.resize(row_idxs.size() + (range.second - range.first + 1ll));
      for(auto i=range.first;i<=range.second;++i,++j)
        row_idxs[j] = i;
    }
    query_config.set_rows_to_query(row_idxs);
  }
  //Attributes
  query_config.set_attributes_to_query(m_attributes);
}

//Loader config functions
JSONLoaderConfig::JSONLoaderConfig(
  bool vid_mapper_file_required) : JSONConfigBase()
{
  m_standalone_converter_process = false;
  m_treat_deletions_as_intervals = false;
  m_produce_combined_vcf = false;
  m_produce_tiledb_array = false;
  m_compress_tiledb_array = true;
  m_disable_synced_writes = false;
  m_delete_and_create_tiledb_array = false;
  m_row_based_partitioning = false;
  //Flag that controls whether the VCF indexes should be discarded to reduce memory consumption
  m_discard_vcf_index = true;
  m_num_entries_in_circular_buffer = 1;
  m_num_converter_processes = 0;
  m_per_partition_size = 0;
  m_max_size_per_callset = 0;
  //Array domain
  m_max_num_rows_in_array = INT64_MAX;
  //#VCF files to open/process in parallel
  m_num_parallel_vcf_files = 1;
  //do ping-pong buffering
  m_do_ping_pong_buffering = true;
  //Offload VCF output processing to another thread
  m_offload_vcf_output_processing = false;
  //Ignore cells that do not belong to this partition
  m_ignore_cells_not_in_partition = false;
  m_vid_mapping_file = "";
  m_callset_mapping_file = "";
  m_segment_size = 10u*1024u*1024u; //10MiB default
  m_num_cells_per_tile = 1024u;
  m_vid_mapper_file_required = vid_mapper_file_required;
  m_fail_if_updating = false;
  m_tiledb_compression_level = Z_DEFAULT_COMPRESSION;
  m_consolidate_tiledb_array_after_load = false;
  m_discard_missing_GTs = false;
  m_no_mandatory_VCF_fields = false;
}

void JSONLoaderConfig::read_from_file(const std::string& filename, VidMapper* id_mapper, const int rank)
{
  JSONConfigBase::read_from_file(filename, id_mapper, rank);
  //Check for row based partitioning - default column based
  m_row_based_partitioning = m_json.HasMember("row_based_partitioning") && m_json["row_based_partitioning"].IsBool()
    && m_json["row_based_partitioning"].GetBool();
  if(m_row_based_partitioning) //Row based partitioning
  {
    VERIFY_OR_THROW(m_json.HasMember("row_partitions"));
  }
  else //Column partitions - if no row based partitioning (default: column partitioning)
  {
    VERIFY_OR_THROW(m_json.HasMember("column_partitions"));
  }
  //Buffer size per column partition
  VERIFY_OR_THROW(m_json.HasMember("size_per_column_partition"));
  m_per_partition_size = m_json["size_per_column_partition"].GetInt64();
  //Obtain number of converters
  m_num_converter_processes = 0;
  if(m_json.HasMember("num_converter_processes") && !m_row_based_partitioning)
    m_num_converter_processes = m_json["num_converter_processes"].GetInt64();
  //Converter processes run independent of loader when num_converter_processes > 0
  m_standalone_converter_process = (m_num_converter_processes && !m_row_based_partitioning) ? true : false;
  //treat deletions as intervals
  if(m_json.HasMember("treat_deletions_as_intervals"))
    m_treat_deletions_as_intervals = m_json["treat_deletions_as_intervals"].GetBool();
  else
    m_treat_deletions_as_intervals = false;
  //Domain size of the array
  m_max_num_rows_in_array = INT64_MAX;
  if(m_json.HasMember("max_num_rows_in_array"))
    m_max_num_rows_in_array = m_json["max_num_rows_in_array"].GetInt64();
  //Ignore callsets with row idx < specified value
  m_lb_callset_row_idx = 0;
  if(m_json.HasMember("lb_callset_row_idx"))
    m_lb_callset_row_idx = m_json["lb_callset_row_idx"].GetInt64();
  //Ignore callsets with row idx > specified value
  m_ub_callset_row_idx = INT64_MAX-1;
  if(m_json.HasMember("ub_callset_row_idx"))
    m_ub_callset_row_idx = m_json["ub_callset_row_idx"].GetInt64();
  //Produce combined vcf
  m_produce_combined_vcf = false;
  if(m_json.HasMember("produce_combined_vcf") && m_json["produce_combined_vcf"].GetBool())
    m_produce_combined_vcf = true;
  //Produce TileDB array
  m_produce_tiledb_array = false;
  if(m_json.HasMember("produce_tiledb_array") && m_json["produce_tiledb_array"].GetBool())
    m_produce_tiledb_array = true;
  //Compress TileDB array by default or if flag set to true
  m_compress_tiledb_array = (!m_json.HasMember("compress_tiledb_array")
      || (m_json["compress_tiledb_array"].IsBool() && m_json["compress_tiledb_array"].GetBool()));
  //Disable synced writes - default false
  m_disable_synced_writes = (m_json.HasMember("disable_synced_writes") && m_json["disable_synced_writes"].IsBool()
      && m_json["disable_synced_writes"].GetBool());
  //recreate array from scratch
  m_delete_and_create_tiledb_array = (m_json.HasMember("delete_and_create_tiledb_array") && m_json["delete_and_create_tiledb_array"].IsBool()
    && m_json["delete_and_create_tiledb_array"].GetBool());
  //Control whether VCF indexes should be discarded to save memory
  m_discard_vcf_index = true;
  if(m_json.HasMember("discard_vcf_index"))
    m_discard_vcf_index = m_json["discard_vcf_index"].GetBool();
  //#vcf files to process in parallel
  m_num_parallel_vcf_files = 1;
  if(m_json.HasMember("num_parallel_vcf_files"))
    m_num_parallel_vcf_files = m_json["num_parallel_vcf_files"].GetInt();
  //do ping pong buffering
  m_do_ping_pong_buffering = true;
  if(m_json.HasMember("do_ping_pong_buffering"))
    m_do_ping_pong_buffering = m_json["do_ping_pong_buffering"].GetBool();
  //Offload VCF output processing
  m_offload_vcf_output_processing = false;
  if(m_json.HasMember("offload_vcf_output_processing"))
    m_offload_vcf_output_processing = m_do_ping_pong_buffering && m_json["offload_vcf_output_processing"].GetBool();
  //Ignore cells that do not belong to this partition
  if(m_json.HasMember("ignore_cells_not_in_partition") && m_json["ignore_cells_not_in_partition"].IsBool())
    m_ignore_cells_not_in_partition = m_json["ignore_cells_not_in_partition"].GetBool();
  //TileDB array segment size
  if(m_json.HasMember("segment_size") && m_json["segment_size"].IsInt64())
    m_segment_size = m_json["segment_size"].GetInt64();
  //TileDB array #cells/tile
  if(m_json.HasMember("num_cells_per_tile") && m_json["num_cells_per_tile"].IsInt64())
    m_num_cells_per_tile = m_json["num_cells_per_tile"].GetInt64();
  if(m_json.HasMember("tiledb_compression_level") && m_json["tiledb_compression_level"].IsInt()) {
    int val = m_json["tiledb_compression_level"].GetInt();
    if ((val < Z_DEFAULT_COMPRESSION) || (val > Z_BEST_COMPRESSION))
       val = Z_DEFAULT_COMPRESSION;
    m_tiledb_compression_level = val;
  }
  //Must have path to vid_mapping_file
  if (m_vid_mapper_file_required) {
    VERIFY_OR_THROW(m_json.HasMember("vid_mapping_file"));
  }
  //flag that causes the loader to fail if this is an update (rather than a fresh load)
  m_fail_if_updating = false;
  if(m_json.HasMember("fail_if_updating") && m_json["fail_if_updating"].IsBool())
    m_fail_if_updating = m_json["fail_if_updating"].GetBool();
  read_and_initialize_vid_and_callset_mapping_if_available(id_mapper, rank);
  //consolidate TileDB array after load - merges fragments
  m_consolidate_tiledb_array_after_load = false;
  if(m_json.HasMember("consolidate_tiledb_array_after_load") && m_json["consolidate_tiledb_array_after_load"].IsBool())
    m_consolidate_tiledb_array_after_load = m_json["consolidate_tiledb_array_after_load"].GetBool();
  //Discard entries with ./. or .|. as the GT field
  m_discard_missing_GTs = false;
  if(m_json.HasMember("discard_missing_GTs") && m_json["discard_missing_GTs"].IsBool())
    m_discard_missing_GTs = m_json["discard_missing_GTs"].GetBool();
  //The array will NOT contain mandatory VCF fields (ref, alt, qual, filter) 
  //if this flag is enabled
  m_no_mandatory_VCF_fields = false;
  if(m_json.HasMember("no_mandatory_VCF_fields") && m_json["no_mandatory_VCF_fields"].IsBool())
    m_no_mandatory_VCF_fields = m_json["no_mandatory_VCF_fields"].GetBool();
}
   
#ifdef HTSDIR

void JSONVCFAdapterConfig::read_from_file(const std::string& filename,
    VCFAdapter& vcf_adapter,
    VidMapper* id_mapper,
    std::string output_format, const int rank,
    const size_t combined_vcf_records_buffer_size_limit)
{
  JSONConfigBase::read_from_file(filename, id_mapper, rank);
  //VCF header filename
  if(m_json.HasMember("vcf_header_filename"))
  {
    const rapidjson::Value& v = m_json["vcf_header_filename"];
    //vcf_header_filename could be an array, one vcf_header_filename location for every rank
    if(v.IsArray())
    {
      VERIFY_OR_THROW(rank < static_cast<int>(v.Size()));
      VERIFY_OR_THROW(v[rank].IsString());
      m_vcf_header_filename = v[rank].GetString();
    }
    else //vcf_header_filename is simply a string
    {
      VERIFY_OR_THROW(v.IsString());
      m_vcf_header_filename = v.GetString();
    }
  }
  //VCF output filename
  if(m_json.HasMember("vcf_output_filename"))
  {
    const rapidjson::Value& v = m_json["vcf_output_filename"];
    //vcf_output_filename could be an array, one vcf_output_filename location for every rank
    if(v.IsArray())
    {
      VERIFY_OR_THROW(rank < static_cast<int>(v.Size()));
      VERIFY_OR_THROW(v[rank].IsString());
      m_vcf_output_filename = v[rank].GetString();
    }
    else //vcf_output_filename is simply a string
    {
      VERIFY_OR_THROW(v.IsString());
      m_vcf_output_filename = v.GetString();
    }
  }
  else
    m_vcf_output_filename = "-";        //stdout
  //VCF output could also be specified in column partitions
  if(m_json.HasMember("column_partitions"))
  {
    //column_partitions_array is an array of the form [ { "begin" : <value> }, {"begin":<value>} ]
    auto& column_partitions_array = m_json["column_partitions"];
    VERIFY_OR_THROW(column_partitions_array.IsArray());
    if(rank < static_cast<int>(column_partitions_array.Size()))
    {
      // {"begin": x}
      auto& curr_partition_info_dict = column_partitions_array[static_cast<rapidjson::SizeType>(rank)];
      VERIFY_OR_THROW(curr_partition_info_dict.IsObject());
      VERIFY_OR_THROW(curr_partition_info_dict.HasMember(JSON_LOADER_PARTITION_INFO_BEGIN_FIELD_NAME));
      if(curr_partition_info_dict.HasMember("vcf_output_filename"))
      {
        VERIFY_OR_THROW(curr_partition_info_dict["vcf_output_filename"].IsString());
        m_vcf_output_filename = curr_partition_info_dict["vcf_output_filename"].GetString();
      }
    }
  }
  //Reference genome
  VERIFY_OR_THROW(m_json.HasMember("reference_genome"));
  {
    const rapidjson::Value& v = m_json["reference_genome"];
    //reference_genome could be an array, one reference_genome location for every rank
    if(v.IsArray())
    {
      VERIFY_OR_THROW(rank < static_cast<int>(v.Size()));
      VERIFY_OR_THROW(v[rank].IsString());
      m_reference_genome = v[rank].GetString();
    }
    else //reference_genome is simply a string
    {
      VERIFY_OR_THROW(v.IsString());
      m_reference_genome = v.GetString();
    }
  }
  //Output format - if arg not empty
  if(output_format == "" && m_json.HasMember("vcf_output_format"))
    output_format = m_json["vcf_output_format"].GetString(); 
  //Limit on max #alt alleles so that PL fields get re-computed
  if(m_json.HasMember("max_diploid_alt_alleles_that_can_be_genotyped"))
    m_max_diploid_alt_alleles_that_can_be_genotyped = m_json["max_diploid_alt_alleles_that_can_be_genotyped"].GetInt();
  else
    m_max_diploid_alt_alleles_that_can_be_genotyped = MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED;
  //Don't produce the full VCF, but determine sites with high allele count
  if(m_json.HasMember("determine_sites_with_max_alleles"))
    m_determine_sites_with_max_alleles = m_json["determine_sites_with_max_alleles"].GetInt();
  else
    m_determine_sites_with_max_alleles = 0;
  //Buffer size for combined vcf records
  if(combined_vcf_records_buffer_size_limit == 0u)
    if(m_json.HasMember("combined_vcf_records_buffer_size_limit"))
      m_combined_vcf_records_buffer_size_limit = m_json["combined_vcf_records_buffer_size_limit"].GetInt64();
    else
      m_combined_vcf_records_buffer_size_limit = DEFAULT_COMBINED_VCF_RECORDS_BUFFER_SIZE;
  else
    m_combined_vcf_records_buffer_size_limit = combined_vcf_records_buffer_size_limit;
  //Cannot be 0
  m_combined_vcf_records_buffer_size_limit = std::max<size_t>(1ull, m_combined_vcf_records_buffer_size_limit);
  //GATK CombineGVCF does not produce GT field by default - option to produce GT
  auto produce_GT_field = (m_json.HasMember("produce_GT_field") && m_json["produce_GT_field"].GetBool());
  //GATK CombineGVCF does not produce FILTER field by default - option to produce FILTER
  auto produce_FILTER_field = (m_json.HasMember("produce_FILTER_field") && m_json["produce_FILTER_field"].GetBool());
  //index output VCF file
  auto index_output_VCF = (m_json.HasMember("index_output_VCF") && m_json["index_output_VCF"].GetBool());
  //sites-only query - doesn't produce any of the FORMAT fields
  auto sites_only_query = (m_json.HasMember("sites_only_query") && m_json["sites_only_query"].GetBool());
  //when producing GT, use the min PL value GT for spanning deletions
  auto produce_GT_with_min_PL_value_for_spanning_deletions = (m_json.HasMember("produce_GT_with_min_PL_value_for_spanning_deletions")
      && m_json["produce_GT_with_min_PL_value_for_spanning_deletions"].GetBool());
  vcf_adapter.initialize(m_reference_genome, m_vcf_header_filename, m_vcf_output_filename, output_format, m_combined_vcf_records_buffer_size_limit,
      produce_GT_field, index_output_VCF, produce_FILTER_field,
      sites_only_query, produce_GT_with_min_PL_value_for_spanning_deletions);
}

void JSONVCFAdapterQueryConfig::read_from_file(const std::string& filename, VariantQueryConfig& query_config,
        VCFAdapter& vcf_adapter, VidMapper* id_mapper,
        std::string output_format, const int rank, const size_t combined_vcf_records_buffer_size_limit)
{
  JSONBasicQueryConfig::read_from_file(filename, query_config, id_mapper, rank);
  JSONVCFAdapterConfig::read_from_file(filename, vcf_adapter, id_mapper,
      output_format, rank, combined_vcf_records_buffer_size_limit);
  query_config.set_sites_only_query(vcf_adapter.sites_only_query());
}
#endif
