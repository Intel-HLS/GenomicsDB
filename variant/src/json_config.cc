//Enable asserts
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "json_config.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw RunConfigException(#X);

void JSONConfigBase::read_from_file(const std::string& filename)
{
  std::ifstream ifs(filename.c_str());
  VERIFY_OR_THROW(ifs.is_open());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  m_json.Parse(str.c_str());
}

void JSONBasicQueryConfig::read_from_file(const std::string& filename, VariantQueryConfig& query_config, int rank)
{
  JSONConfigBase::read_from_file(filename);
  //Workspace
  VERIFY_OR_THROW(m_json.HasMember("workspace"));
  {
    const rapidjson::Value& workspace = m_json["workspace"];
    //workspace could be an array, one workspace dir for every rank
    if(workspace.IsArray())
    {
      VERIFY_OR_THROW(rank < workspace.Size());
      VERIFY_OR_THROW(workspace[rank].IsString());
      m_workspace = workspace[rank].GetString();
    }
    else //workspace is simply a string
    {
      VERIFY_OR_THROW(workspace.IsString());
      m_workspace = workspace.GetString();
    }
  }
  //Array
  VERIFY_OR_THROW(m_json.HasMember("array"));
  {
    const rapidjson::Value& array_name = m_json["array"];
    //array could be an array, one array dir for every rank
    if(array_name.IsArray())
    {
      VERIFY_OR_THROW(rank < array_name.Size());
      VERIFY_OR_THROW(array_name[rank].IsString());
      m_array_name = array_name[rank].GetString();
    }
    else //array is simply a string
    {
      VERIFY_OR_THROW(array_name.IsString());
      m_array_name = array_name.GetString();
    }
  }
  //Query columns
  //Example:  [ [ [0,5], 45 ], [ 76, 87 ] ]
  //This means that rank 0 will have 2 query intervals: [0-5] and [45-45] and rank 1 will have
  //2 intervals [76-76] and [87-87]
  //But you could have a single innermost list - with this option all ranks will query the same list 
  VERIFY_OR_THROW(m_json.HasMember("query_column_ranges") || m_json.HasMember("scan_full"));
  if(m_json.HasMember("query_column_ranges"))
  {
    const rapidjson::Value& q1 = m_json["query_column_ranges"];
    VERIFY_OR_THROW(q1.IsArray());
    //Single element list or rank < size
    VERIFY_OR_THROW(q1.Size() == 1 || rank < q1.Size());
    const rapidjson::Value& q2 = (q1.Size() == 1) ? q1[0u] : q1[rank];
    VERIFY_OR_THROW(q2.IsArray());
    for(rapidjson::SizeType i=0;i<q2.Size();++i)
    {
      const rapidjson::Value& q3 = q2[i];
      //q3 is list of 2 elements to represent query interval
      if(q3.IsArray())
      {
        VERIFY_OR_THROW(q3.Size() == 2);
        VERIFY_OR_THROW(q3[0u].IsInt64());
        VERIFY_OR_THROW(q3[1u].IsInt64());
        query_config.add_column_interval_to_query(q3[0u].GetInt64(), q3[1u].GetInt64());
      }
      else //single position
      {
        VERIFY_OR_THROW(q3.IsInt64());
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
    VERIFY_OR_THROW(q1.IsArray());
    //Single element list or rank < size
    VERIFY_OR_THROW(q1.Size() == 1 || rank < q1.Size());
    const rapidjson::Value& q2 = (q1.Size() == 1) ? q1[0u] : q1[rank];
    VERIFY_OR_THROW(q2.IsArray());
    std::vector<int64_t> row_idxs;
    for(rapidjson::SizeType i=0;i<q2.Size();++i)
    {
      const rapidjson::Value& q3 = q2[i];
      //q3 is list of 2 elements to represent query row interval
      if(q3.IsArray())
      {
        VERIFY_OR_THROW(q3.Size() == 2);
        VERIFY_OR_THROW(q3[0u].IsInt64());
        VERIFY_OR_THROW(q3[1u].IsInt64());
        for(auto i=q3[0u].GetInt64();i<=q3[1u].GetInt64();++i)
          row_idxs.push_back(i);
      }
      else //single position
      {
        VERIFY_OR_THROW(q3.IsInt64());
        row_idxs.push_back(q3.GetInt64());
      }
    }
    query_config.set_rows_to_query(row_idxs);
  }
  VERIFY_OR_THROW(m_json.HasMember("query_attributes"));
  {
    const rapidjson::Value& q1 = m_json["query_attributes"];
    VERIFY_OR_THROW(q1.IsArray());
    std::vector<std::string> attributes(q1.Size());
    for(rapidjson::SizeType i=0;i<q1.Size();++i)
    {
      const rapidjson::Value& q2 = q1[i];
      VERIFY_OR_THROW(q2.IsString());
      attributes[i] = std::move(std::string(q2.GetString()));
    }
    query_config.set_attributes_to_query(attributes);
  }
}
   
#ifdef HTSDIR

void JSONVCFAdapterConfig::read_from_file(const std::string& filename,
    VCFAdapter& vcf_adapter, std::string output_format, int rank)
{
  JSONConfigBase::read_from_file(filename);
  //VCF header filename
  VERIFY_OR_THROW(m_json.HasMember("vcf_header_filename"));
  {
    const rapidjson::Value& v = m_json["vcf_header_filename"];
    //vcf_header_filename could be an array, one vcf_header_filename location for every rank
    if(v.IsArray())
    {
      VERIFY_OR_THROW(rank < v.Size());
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
      VERIFY_OR_THROW(rank < v.Size());
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
  //Reference genome
  VERIFY_OR_THROW(m_json.HasMember("reference_genome"));
  {
    const rapidjson::Value& v = m_json["reference_genome"];
    //reference_genome could be an array, one reference_genome location for every rank
    if(v.IsArray())
    {
      VERIFY_OR_THROW(rank < v.Size());
      VERIFY_OR_THROW(v[rank].IsString());
      m_reference_genome = v[rank].GetString();
    }
    else //reference_genome is simply a string
    {
      VERIFY_OR_THROW(v.IsString());
      m_reference_genome = v.GetString();
    }
  }
  vcf_adapter.initialize(m_reference_genome, m_vcf_header_filename, m_vcf_output_filename, output_format);
}

void JSONVCFAdapterQueryConfig::read_from_file(const std::string& filename, VariantQueryConfig& query_config,
        VCFAdapter& vcf_adapter, FileBasedVidMapper& id_mapper,
        std::string output_format, int rank)
{
  JSONBasicQueryConfig::read_from_file(filename, query_config, rank);
  JSONVCFAdapterConfig::read_from_file(filename, vcf_adapter, output_format, rank);
  //Over-ride callset mapping file in top-level config if necessary
  std::string callset_mapping_file="";
  if(JSONBasicQueryConfig::m_json.HasMember("callset_mapping_file") &&
      JSONBasicQueryConfig::m_json["callset_mapping_file"].IsString())
    callset_mapping_file = JSONBasicQueryConfig::m_json["callset_mapping_file"].GetString();
  //contig and callset id mapping
  VERIFY_OR_THROW(JSONBasicQueryConfig::m_json.HasMember("vid_mapping_file"));
  {
    const rapidjson::Value& v = JSONBasicQueryConfig::m_json["vid_mapping_file"];
    //Could be array - one for each process
    if(v.IsArray())
    {
      VERIFY_OR_THROW(rank < v.Size());
      VERIFY_OR_THROW(v[rank].IsString());
      id_mapper = std::move(FileBasedVidMapper(v[rank].GetString(), callset_mapping_file));
    }
    else //or single string for all processes
    {
      VERIFY_OR_THROW(v.IsString());
      id_mapper = std::move(FileBasedVidMapper(v.GetString(), callset_mapping_file));
    }
  }
}
#endif
