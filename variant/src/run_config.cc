#include "run_config.h"

//Enable asserts
#ifdef NDEBUG
#undef NDEBUG
#endif

RunConfig g_run_config;

void RunConfig::read_from_file(const std::string& filename, VariantQueryConfig& query_config, int rank)
{
  std::ifstream ifs(filename.c_str());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  m_json.Parse(str.c_str());
  //Workspace
  assert(m_json.HasMember("workspace"));
  {
    const rapidjson::Value& workspace = m_json["workspace"];
    //workspace could be an array, one workspace dir for every rank
    if(workspace.IsArray())
    {
      assert(rank < workspace.Size());
      assert(workspace[rank].IsString());
      m_workspace = workspace[rank].GetString();
    }
    else //workspace is simply a string
    {
      assert(workspace.IsString());
      m_workspace = workspace.GetString();
    }
  }
  //Array
  assert(m_json.HasMember("array"));
  {
    const rapidjson::Value& array_name = m_json["array"];
    //array could be an array, one array dir for every rank
    if(array_name.IsArray())
    {
      assert(rank < array_name.Size());
      assert(array_name[rank].IsString());
      m_array_name = array_name[rank].GetString();
    }
    else //array is simply a string
    {
      assert(array_name.IsString());
      m_array_name = array_name.GetString();
    }
  }
  //Query columns
  //Example:  [ [ [0,5], 45 ], [ 76, 87 ] ]
  //This means that rank 0 will have 2 query intervals: [0-5] and [45-45] and rank 1 will have
  //2 intervals [76-76] and [87-87]
  //But you could have a single innermost list - with this option all ranks will query the same list 
  assert(m_json.HasMember("query_column_ranges"));
  {
    const rapidjson::Value& q1 = m_json["query_column_ranges"];
    assert(q1.IsArray());
    //Single element list or rank < size
    assert(q1.Size() == 1 || rank < q1.Size());
    const rapidjson::Value& q2 = (q1.Size() == 1) ? q1[0u] : q1[rank];
    assert(q2.IsArray());
    for(rapidjson::SizeType i=0;i<q2.Size();++i)
    {
      const rapidjson::Value& q3 = q2[i];
      //q3 is list of 2 elements to represent query interval
      if(q3.IsArray())
      {
        assert(q3.Size() == 2);
        assert(q3[0u].IsInt64());
        assert(q3[1u].IsInt64());
        query_config.add_column_interval_to_query(q3[0u].GetInt64(), q3[1u].GetInt64());
      }
      else //single position
      {
        assert(q3.IsInt64());
        query_config.add_column_interval_to_query(q3.GetInt64(), q3.GetInt64());
      }
    }
  }
  //Query rows
  //Example:  [ [ [0,5], 45 ], [ 76, 87 ] ]
  //This means that rank 0 will query rows: [0-5] and [45-45] and rank 1 will have
  //2 intervals [76-76] and [87-87]
  //But you could have a single innermost list - with this option all ranks will query the same list 
  if(m_json.HasMember("query_row_ranges"))
  {
    const rapidjson::Value& q1 = m_json["query_row_ranges"];
    assert(q1.IsArray());
    //Single element list or rank < size
    assert(q1.Size() == 1 || rank < q1.Size());
    const rapidjson::Value& q2 = (q1.Size() == 1) ? q1[0u] : q1[rank];
    assert(q2.IsArray());
    std::vector<int64_t> row_idxs;
    for(rapidjson::SizeType i=0;i<q2.Size();++i)
    {
      const rapidjson::Value& q3 = q2[i];
      //q3 is list of 2 elements to represent query row interval
      if(q3.IsArray())
      {
        assert(q3.Size() == 2);
        assert(q3[0u].IsInt64());
        assert(q3[1u].IsInt64());
        for(auto i=q3[0u].GetInt64();i<=q3[1u].GetInt64();++i)
          row_idxs.push_back(i);
      }
      else //single position
      {
        assert(q3.IsInt64());
        row_idxs.push_back(q3.GetInt64());
      }
    }
    query_config.set_rows_to_query(row_idxs);
  }
  assert(m_json.HasMember("query_attributes"));
  {
    const rapidjson::Value& q1 = m_json["query_attributes"];
    assert(q1.IsArray());
    std::vector<std::string> attributes(q1.Size());
    for(rapidjson::SizeType i=0;i<q1.Size();++i)
    {
      const rapidjson::Value& q2 = q1[i];
      assert(q2.IsString());
      attributes[i] = std::move(std::string(q2.GetString()));
    }
    query_config.set_attributes_to_query(attributes);
  }
}
