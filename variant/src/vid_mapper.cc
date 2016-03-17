#include "vid_mapper.h"
#include "rapidjson/document.h"
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"
#include "c_api.h"

std::unordered_map<std::string, int> VidMapper::m_length_descriptor_string_to_int = std::unordered_map<std::string, int>({
    {"BCF_VL_FIXED", BCF_VL_FIXED},
    {"BCF_VL_A", BCF_VL_A},
    {"A", BCF_VL_A},
    {"BCF_VL_R", BCF_VL_R},
    {"R", BCF_VL_R},
    {"BCF_VL_G", BCF_VL_G},
    {"G", BCF_VL_G},
    {"BCF_VL_P", BCF_VL_P},
    {"P", BCF_VL_P},
    {"BCF_VL_VAR", BCF_VL_VAR},
    {"VAR", BCF_VL_VAR}
    });

std::unordered_map<std::string, std::type_index> VidMapper::m_typename_string_to_type_index =
  std::unordered_map<std::string, std::type_index>({
      {"int", std::type_index(typeid(int))},
      {"float", std::type_index(typeid(float))},
      {"char", std::type_index(typeid(char))}
      });

std::unordered_map<std::string, int> VidMapper::m_typename_string_to_bcf_ht_type =
  std::unordered_map<std::string, int>({
      {"int", BCF_HT_INT},
      {"float", BCF_HT_REAL},
      {"char", BCF_HT_STR}
      });

#define VERIFY_OR_THROW(X) if(!(X)) throw VidMapperException(#X);

void VidMapper::clear()
{
  m_callset_name_to_row_idx.clear();
  m_row_idx_to_info.clear();
  m_filename_to_idx.clear();
  m_file_idx_to_info.clear();
  m_contig_name_to_idx.clear();
  m_contig_idx_to_info.clear();
  m_contig_begin_2_idx.clear();
  m_contig_end_2_idx.clear();
  m_field_name_to_idx.clear();
  m_field_idx_to_info.clear();
  m_owner_idx_to_file_idx_vec.clear();
}

bool VidMapper::get_contig_location(int64_t query_position, std::string& contig_name, int64_t& contig_position) const
{
  int idx = -1;
  std::pair<int64_t, int> query_pair;
  query_pair.first = query_position;
  query_pair.second = 0;
  //find contig with offset >= query_position
  auto iter = std::lower_bound(m_contig_begin_2_idx.begin(), m_contig_begin_2_idx.end(), query_pair, contig_offset_idx_pair_cmp);
  if(iter == m_contig_begin_2_idx.end())        //no such contig exists, hence, get last contig in sorted order
  {
    assert(m_contig_begin_2_idx.size() > 0u);
    idx = m_contig_begin_2_idx[m_contig_begin_2_idx.size()-1u].second;
  }
  else
  {
    if((*iter).first == query_position)       //query_position == contig offset at iter, found idx
      idx = (*iter).second;
    else                                //query_position < contig_offset at iter, get idx at previous element
    {
      //if iter == begin(), query position is less than 1st contig offset, return invalid
      if(iter == m_contig_begin_2_idx.begin()) 
        return false;
      auto vector_idx = iter - m_contig_begin_2_idx.begin();
      idx = m_contig_begin_2_idx[vector_idx-1].second;
    }
  }
  if(idx < 0)
    return false;
  assert(static_cast<size_t>(idx) < m_contig_idx_to_info.size());
  //query_position is within the contig
  auto contig_offset = m_contig_idx_to_info[idx].m_tiledb_column_offset;
  auto contig_length = m_contig_idx_to_info[idx].m_length;
  if((query_position >= contig_offset) && (query_position < contig_offset+contig_length))
  {
    contig_name = m_contig_idx_to_info[idx].m_name;
    contig_position = query_position - contig_offset;
    return true;
  }
  return false;
}

bool VidMapper::get_next_contig_location(int64_t query_position, std::string& next_contig_name, int64_t& next_contig_offset) const
{
  int idx = -1;
  std::pair<int64_t, int> query_pair;
  query_pair.first = query_position;
  query_pair.second = 0;
  //find contig with offset > query_position
  auto iter = std::upper_bound(m_contig_begin_2_idx.begin(), m_contig_begin_2_idx.end(), query_pair, contig_offset_idx_pair_cmp);
  if(iter == m_contig_begin_2_idx.end())        //no such contig exists, hence, set large upper bound
  {
    next_contig_name = "";
    next_contig_offset = INT64_MAX;
    return false;
  }
  else
  {
    idx = (*iter).second;
    assert(idx >=0 && static_cast<size_t>(idx) < m_contig_idx_to_info.size());
    next_contig_name = m_contig_idx_to_info[idx].m_name;
    next_contig_offset = m_contig_idx_to_info[idx].m_tiledb_column_offset;
    assert(next_contig_offset > query_position);
    return true;
  }
}

bool VidMapper::get_tiledb_position(int64_t& position, const std::string& contig_name, const int64_t contig_position) const
{
  auto iter = m_contig_name_to_idx.find(contig_name);
  if(iter == m_contig_name_to_idx.end())
    return false;
  auto idx = (*iter).second;
  assert(static_cast<size_t>(idx) < m_contig_idx_to_info.size());
  if(contig_position >= m_contig_idx_to_info[idx].m_length)
    return false;
  position = m_contig_idx_to_info[idx].m_tiledb_column_offset + contig_position;
  return true;
}
    
bool VidMapper::get_callset_name(const int64_t row_idx, std::string& callset_name) const
{
  if(static_cast<size_t>(row_idx) >= m_row_idx_to_info.size())
    return false;
  callset_name = m_row_idx_to_info[row_idx].m_name;
  return true;
}
    
bool VidMapper::get_tiledb_row_idx(int64_t& row_idx, const std::string& callset_name) const
{
  auto iter = m_callset_name_to_row_idx.find(callset_name);
  if(iter == m_callset_name_to_row_idx.end())
    return false;
  row_idx = (*iter).second;
  return true;
}

void VidMapper::build_vcf_fields_vectors(std::vector<std::vector<std::string>>& vcf_fields) const
{
  vcf_fields.clear();
  vcf_fields.resize(3u);        //FILTER,INFO,FORMAT
  for(const auto& field_info : m_field_idx_to_info)
  {
    if(field_info.m_is_vcf_INFO_field)
      vcf_fields[BCF_HL_INFO].push_back(field_info.m_name);
    if(field_info.m_is_vcf_FORMAT_field)
      vcf_fields[BCF_HL_FMT].push_back(field_info.m_name);
  }
}

void VidMapper::build_tiledb_array_schema(VariantArraySchema*& array_schema, const std::string array_name,
    const bool row_based_partitioning, const RowRange& row_range, const bool compress_fields)
  const
{
  auto dim_names = std::vector<std::string>({"samples", "position"});
  auto dim_domains = row_based_partitioning ?
    std::vector<std::pair<int64_t, int64_t>>({ {row_range.first, row_range.second}, {0, INT64_MAX}})
    : std::vector<std::pair<int64_t,int64_t>>({ {0, get_num_callsets()-1}, {0, INT64_MAX } });
  std::vector<std::string> attribute_names;
  std::vector<std::type_index> types;
  std::vector<int> num_vals;
  //END
  attribute_names.push_back("END");
  types.push_back(std::type_index(typeid(int64_t)));
  num_vals.push_back(1);
  //REF
  attribute_names.push_back("REF");
  types.push_back(std::type_index(typeid(char)));
  num_vals.push_back(TILEDB_VAR_NUM);
  //ALT
  attribute_names.push_back("ALT");
  types.push_back(std::type_index(typeid(char)));
  num_vals.push_back(TILEDB_VAR_NUM);
  //QUAL
  attribute_names.push_back("QUAL");
  types.push_back(std::type_index(typeid(float)));
  num_vals.push_back(1);
  //FILTER
  attribute_names.push_back("FILTER");
  types.push_back(std::type_index(typeid(int)));
  num_vals.push_back(TILEDB_VAR_NUM);
  //INFO fields
  for(const auto& field_info : m_field_idx_to_info)
  {
    if(field_info.m_name == "END")      //skip END field
      continue;
    if(field_info.m_is_vcf_INFO_field)
    {
      attribute_names.push_back(field_info.m_name);
      types.push_back(field_info.m_type_index);
      num_vals.push_back(field_info.m_length_descriptor == BCF_VL_FIXED ? field_info.m_num_elements : TILEDB_VAR_NUM);
    }
  }
  //FORMAT fields
  for(const auto& field_info : m_field_idx_to_info)
  {
    if(field_info.m_name == "END")      //skip END field
      continue;
    if(field_info.m_is_vcf_FORMAT_field)
    {
      if(field_info.m_is_vcf_INFO_field)        //Also an INFO field of the same name - add suffix
        attribute_names.push_back(field_info.m_name+"_FORMAT");
      else
        attribute_names.push_back(field_info.m_name);
      types.push_back(field_info.m_type_index);
      num_vals.push_back(field_info.m_length_descriptor == BCF_VL_FIXED ? field_info.m_num_elements : TILEDB_VAR_NUM);
    }
  }
  //COORDS
  types.push_back(std::type_index(typeid(int64_t)));
  //For compression
  //no compression - empty vector
  std::vector<int> compression;
  if(compress_fields)
    for(auto i=0u;i<types.size();++i)   //types contains entry for coords also
      compression.push_back(TILEDB_GZIP);
  else
    for(auto i=0u;i<types.size();++i)   //types contains entry for coords also
      compression.push_back(TILEDB_NO_COMPRESSION);
  array_schema = new VariantArraySchema(array_name, attribute_names, dim_names, dim_domains, types, num_vals, compression);
}

void VidMapper::build_file_partitioning(const int partition_idx, const RowRange row_partition)
{
  if(static_cast<size_t>(partition_idx) >= m_owner_idx_to_file_idx_vec.size())
    m_owner_idx_to_file_idx_vec.resize(partition_idx+1);
  m_owner_idx_to_file_idx_vec[partition_idx].clear();
  auto max_row_idx = std::min<int64_t>(row_partition.second, get_num_callsets()-1);
  std::unordered_set<int64_t> files_set;
  for(auto row_idx=row_partition.first;row_idx<=max_row_idx;++row_idx)
  {
    VERIFY_OR_THROW(row_idx >= 0 && static_cast<size_t>(row_idx) < m_row_idx_to_info.size());
    auto file_idx = m_row_idx_to_info[row_idx].m_file_idx;
    VERIFY_OR_THROW(file_idx >= 0 && static_cast<size_t>(file_idx) < m_file_idx_to_info.size());
    if(files_set.find(file_idx) == files_set.end())
    {
      auto& file_info = m_file_idx_to_info[file_idx];
      file_info.m_owner_idx = partition_idx;
      m_owner_idx_to_file_idx_vec[partition_idx].push_back(file_idx);
      auto& local_row_idx_pairs_vec = file_info.m_local_tiledb_row_idx_pairs;
      auto last_valid_idx = -1ll;
      //Check if any callset within this file is outside the row partition
      for(auto i=0ull;i<local_row_idx_pairs_vec.size();++i)
      {
        //within range
        if(local_row_idx_pairs_vec[i].second >= row_partition.first && local_row_idx_pairs_vec[i].second <= row_partition.second)
        {
          ++last_valid_idx;
          if(i > static_cast<uint64_t>(last_valid_idx))
            std::swap(local_row_idx_pairs_vec[last_valid_idx], local_row_idx_pairs_vec[i]);
        }
      }
      local_row_idx_pairs_vec.resize(last_valid_idx+1);
      files_set.insert(file_idx);
    }
  }
  sort_and_assign_local_file_idxs_for_partition(partition_idx);
}

void VidMapper::sort_and_assign_local_file_idxs_for_partition(const int owner_idx)
{
  //Sort files by global_file_idx
  std::sort(m_owner_idx_to_file_idx_vec[owner_idx].begin(), m_owner_idx_to_file_idx_vec[owner_idx].end());
  //Set local idx in file idx to info struct
  for(auto i=0u;i<m_owner_idx_to_file_idx_vec[owner_idx].size();++i)
  {
    auto global_file_idx = m_owner_idx_to_file_idx_vec[owner_idx][i];
    assert(static_cast<size_t>(global_file_idx) < m_file_idx_to_info.size());
    m_file_idx_to_info[global_file_idx].m_local_file_idx = i;
  }
}

void VidMapper::verify_file_partitioning() const
{
  auto unassigned_files = false;
  for(auto file_idx=0ull;file_idx<m_file_idx_to_info.size();++file_idx)
  {
    auto& file_info = m_file_idx_to_info[file_idx];
    if(file_info.m_owner_idx < 0)
    {
      std::cerr << "File " << file_info.m_name << " is not assigned to any partition\n";
      unassigned_files = true;
    }
  }
  if(unassigned_files)
    throw FileBasedVidMapperException("Found files that are not assigned to any partition");
}

//FileBasedVidMapper code
#ifdef VERIFY_OR_THROW
#undef VERIFY_OR_THROW
#endif

#define VERIFY_OR_THROW(X) if(!(X)) throw FileBasedVidMapperException(#X);

FileBasedVidMapper::FileBasedVidMapper(const std::string& filename, const std::string& callset_mapping_file,
    const int64_t limit_callset_row_idx, const bool callsets_file_required)
  : VidMapper()
{
  VERIFY_OR_THROW(filename.length() && "Vid mapping file unspecified");
  rapidjson::Document json_doc;
  std::ifstream ifs(filename.c_str());
  if(!ifs.is_open())
    throw FileBasedVidMapperException((std::string("Could not open vid mapping file \"")+filename+"\"").c_str());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  json_doc.Parse(str.c_str());
  m_limit_callset_row_idx = limit_callset_row_idx;
  //Callset info parsing
  if(callset_mapping_file != "")
    parse_callsets_file(callset_mapping_file);
  else
  {
    VERIFY_OR_THROW(!callsets_file_required || (json_doc.HasMember("callset_mapping_file") && json_doc["callset_mapping_file"].IsString()));
    if(callsets_file_required || json_doc.HasMember("callset_mapping_file"))
      parse_callsets_file(json_doc["callset_mapping_file"].GetString());
  }
  //Contig info parsing
  VERIFY_OR_THROW(json_doc.HasMember("contigs"));
  {
    const rapidjson::Value& contigs_dict = json_doc["contigs"];
    //contigs is a dictionary of name:info key-value pairs
    VERIFY_OR_THROW(contigs_dict.IsObject());
    auto num_contigs = contigs_dict.MemberCount();
    m_contig_idx_to_info.resize(num_contigs);
    m_contig_begin_2_idx.resize(num_contigs);
    m_contig_end_2_idx.resize(num_contigs);
    std::string contig_name;
    std::string filename;
    auto contig_idx=0;
    for(auto b=contigs_dict.MemberBegin(), e=contigs_dict.MemberEnd();b!=e;++b,++contig_idx)
    {
      const auto& curr_obj = *b;
      contig_name = curr_obj.name.GetString();
      const auto& contig_info_dict = curr_obj.value;
      VERIFY_OR_THROW(contig_info_dict.IsObject());    //must be dict
      VERIFY_OR_THROW(contig_info_dict.HasMember("tiledb_column_offset") && contig_info_dict["tiledb_column_offset"].IsInt64());
      auto tiledb_column_offset = contig_info_dict["tiledb_column_offset"].GetInt64();
      VERIFY_OR_THROW(contig_info_dict.HasMember("length") && contig_info_dict["length"].IsInt64());
      auto length = contig_info_dict["length"].GetInt64();
      VERIFY_OR_THROW(static_cast<size_t>(contig_idx) < static_cast<size_t>(num_contigs));
      m_contig_name_to_idx[contig_name] = contig_idx;
      m_contig_idx_to_info[contig_idx].set_info(contig_idx, contig_name, length, tiledb_column_offset);
      m_contig_begin_2_idx[contig_idx].first = tiledb_column_offset; 
      m_contig_begin_2_idx[contig_idx].second = contig_idx; 
      m_contig_end_2_idx[contig_idx].first = tiledb_column_offset + length - 1; //inclusive
      m_contig_end_2_idx[contig_idx].second = contig_idx; 
    }
    std::sort(m_contig_begin_2_idx.begin(), m_contig_begin_2_idx.end(), contig_offset_idx_pair_cmp);
    std::sort(m_contig_end_2_idx.begin(), m_contig_end_2_idx.end(), contig_offset_idx_pair_cmp);
  }
  //Field info parsing
  VERIFY_OR_THROW(json_doc.HasMember("fields"));
  {
    const rapidjson::Value& fields_dict = json_doc["fields"];
    //fields is a dictionary: name: dict
    VERIFY_OR_THROW(fields_dict.IsObject());
    auto num_fields = fields_dict.MemberCount();
    m_field_idx_to_info.resize(num_fields);
    std::string field_name;
    auto field_idx = 0;
    for(auto b=fields_dict.MemberBegin(), e=fields_dict.MemberEnd();b!=e;++b,++field_idx)
    {
      const auto& curr_obj = *b;
      field_name = curr_obj.name.GetString();
      m_field_name_to_idx[field_name] = field_idx;
      m_field_idx_to_info[field_idx].set_info(field_name, field_idx);
      const auto& field_info_dict = curr_obj.value;
      VERIFY_OR_THROW(field_info_dict.IsObject());
      //Field type - int, char etc
      VERIFY_OR_THROW(field_info_dict.HasMember("type") && field_info_dict["type"].IsString());
      {
        auto iter = VidMapper::m_typename_string_to_type_index.find(field_info_dict["type"].GetString());
        VERIFY_OR_THROW(iter != VidMapper::m_typename_string_to_type_index.end() && "Unhandled field type");
        m_field_idx_to_info[field_idx].m_type_index = (*iter).second;
      }
      {
        auto iter = VidMapper::m_typename_string_to_bcf_ht_type.find(field_info_dict["type"].GetString());
        VERIFY_OR_THROW(iter != VidMapper::m_typename_string_to_bcf_ht_type.end() && "Unhandled field type");
        m_field_idx_to_info[field_idx].m_bcf_ht_type = (*iter).second;
      }
      if(field_info_dict.HasMember("vcf_field_class"))
      {
        //Array which specifies whether field if INFO, FORMAT, FILTER etc
        const auto& vcf_field_class_array = field_info_dict["vcf_field_class"];
        VERIFY_OR_THROW(vcf_field_class_array.IsArray());
        for(rapidjson::SizeType i=0;i<vcf_field_class_array.Size();++i)
        {
          auto class_name = std::move(std::string(vcf_field_class_array[i].GetString()));
          if(class_name == "INFO")
            m_field_idx_to_info[field_idx].m_is_vcf_INFO_field = true;
          else
            if(class_name == "FORMAT")
              m_field_idx_to_info[field_idx].m_is_vcf_FORMAT_field = true;
        }
      }
      if(field_info_dict.HasMember("length"))
      {
        if(field_info_dict["length"].IsInt64())
          m_field_idx_to_info[field_idx].m_num_elements = field_info_dict["length"].GetInt64();
        else
        {
          VERIFY_OR_THROW(field_info_dict["length"].IsString());
          auto iter = VidMapper::m_length_descriptor_string_to_int.find(field_info_dict["length"].GetString());
          if(iter == VidMapper::m_length_descriptor_string_to_int.end())
            m_field_idx_to_info[field_idx].m_length_descriptor = BCF_VL_VAR;
          else
            m_field_idx_to_info[field_idx].m_length_descriptor = (*iter).second;
        }
      }
    }
  } 
  m_is_initialized = true;
}

void FileBasedVidMapper::parse_callsets_file(const std::string& filename)
{
  VERIFY_OR_THROW(filename.length() && "Vid mapping file unspecified");
  rapidjson::Document json_doc;
  std::ifstream ifs(filename.c_str());
  if(!ifs.is_open())
    throw FileBasedVidMapperException((std::string("Could not open callsets file \"")+filename+"\"").c_str());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  json_doc.Parse(str.c_str());
  uint64_t num_files = 0ull;
  //Callset info parsing
  VERIFY_OR_THROW(json_doc.HasMember("callsets"));
  {
    const rapidjson::Value& callsets_dict = json_doc["callsets"];
    //callsets is a dictionary of name:info key-value pairs
    VERIFY_OR_THROW(callsets_dict.IsObject());
    auto num_callsets = (m_limit_callset_row_idx == INT64_MAX) ? callsets_dict.MemberCount() :
      std::min<int64_t>(callsets_dict.MemberCount(), m_limit_callset_row_idx+1);
    m_row_idx_to_info.resize(num_callsets);
    std::string callset_name;
    for(auto b=callsets_dict.MemberBegin(), e=callsets_dict.MemberEnd();b!=e;++b)
    {
      const auto& curr_obj = *b;
      callset_name = curr_obj.name.GetString();
      const auto& callset_info_dict = curr_obj.value;
      VERIFY_OR_THROW(callset_info_dict.IsObject());    //must be dict
      VERIFY_OR_THROW(callset_info_dict.HasMember("row_idx"));
      int64_t row_idx = callset_info_dict["row_idx"].GetInt64();
      if(row_idx > m_limit_callset_row_idx)
        continue;
      int64_t file_idx = -1;
      if(callset_info_dict.HasMember("filename"))
      {
        std::string filename = std::move(callset_info_dict["filename"].GetString());
        auto iter = m_filename_to_idx.find(filename);
        if(iter == m_filename_to_idx.end())
        {
          iter = m_filename_to_idx.insert(std::make_pair(filename, num_files++)).first;
          file_idx = num_files-1;
          if(num_files > m_file_idx_to_info.size())
            m_file_idx_to_info.resize(2*num_files+1);
          m_file_idx_to_info[file_idx].set_info(file_idx, filename);
        }
        else
          file_idx = (*iter).second;
        //local callset idx
        auto idx_in_file = 0ll;
        if(callset_info_dict.HasMember("idx_in_file"))
          idx_in_file = callset_info_dict["idx_in_file"].GetInt64();
        assert(file_idx < static_cast<int64_t>(m_file_idx_to_info.size()));
        m_file_idx_to_info[file_idx].add_local_tiledb_row_idx_pair(idx_in_file, row_idx);
      }
      m_callset_name_to_row_idx[callset_name] = row_idx;
      VERIFY_OR_THROW(static_cast<size_t>(row_idx) < m_row_idx_to_info.size());
      m_row_idx_to_info[row_idx].set_info(row_idx, callset_name, file_idx);
    }
  }
  m_file_idx_to_info.resize(num_files);
  //File partitioning info
  if(json_doc.HasMember("file_division"))
  {
    const auto& file_division_array = json_doc["file_division"];
    //is an array
    VERIFY_OR_THROW(file_division_array.IsArray());
    m_owner_idx_to_file_idx_vec.resize(file_division_array.Size());
    for(rapidjson::SizeType owner_idx=0;owner_idx<file_division_array.Size();++owner_idx)
    {
      auto& files_array = file_division_array[owner_idx];
      VERIFY_OR_THROW(files_array.IsArray());
      std::unordered_set<int64_t> files_set;
      for(rapidjson::SizeType i=0;i<files_array.Size();++i)
      {
        int64_t global_file_idx;
        auto found = get_global_file_idx(files_array[i].GetString(), global_file_idx);
        if(found && files_set.find(global_file_idx) == files_set.end())
        {
          m_file_idx_to_info[global_file_idx].m_owner_idx = owner_idx;
          m_owner_idx_to_file_idx_vec[owner_idx].push_back(global_file_idx);
          files_set.insert(global_file_idx);
        }
      }
      sort_and_assign_local_file_idxs_for_partition(owner_idx);
    }
  }
}
