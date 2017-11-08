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

#include <libgen.h>
#include "vid_mapper.h"
#include "tiledb.h"
#include "json_config.h"
#include "known_field_info.h"
#include "vcf.h"

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
    {"VAR", BCF_VL_VAR},
    {"PP", BCF_VL_Phased_Ploidy},
    {"Phased_Ploidy", BCF_VL_Phased_Ploidy},
    {"phased_ploidy", BCF_VL_Phased_Ploidy},
    {"PHASED_PLOIDY", BCF_VL_Phased_Ploidy}
    });

std::unordered_map<std::string, std::type_index> VidMapper::m_typename_string_to_type_index =
  std::unordered_map<std::string, std::type_index>({
      {"int", std::type_index(typeid(int))},
      {"Int", std::type_index(typeid(int))},
      {"integer", std::type_index(typeid(int))},
      {"Integer", std::type_index(typeid(int))},
      {"float", std::type_index(typeid(float))},
      {"Float", std::type_index(typeid(float))}, 
      {"bool", std::type_index(typeid(char))},
      {"Bool", std::type_index(typeid(char))},
      {"boolean", std::type_index(typeid(char))},
      {"Boolean", std::type_index(typeid(char))},
      {"flag", std::type_index(typeid(char))},
      {"Flag", std::type_index(typeid(char))},
      {"string", std::type_index(typeid(char))},
      {"String", std::type_index(typeid(char))},
      {"char", std::type_index(typeid(char))},
      {"Char", std::type_index(typeid(char))}
      });

std::unordered_map<std::string, int> VidMapper::m_typename_string_to_bcf_ht_type =
  std::unordered_map<std::string, int>({
      {"int", BCF_HT_INT},
      {"Int", BCF_HT_INT},
      {"integer", BCF_HT_INT},
      {"Integer", BCF_HT_INT},
      {"float", BCF_HT_REAL},
      {"Float", BCF_HT_REAL},
      {"bool", BCF_HT_FLAG},
      {"Bool", BCF_HT_FLAG},
      {"boolean", BCF_HT_FLAG},
      {"Boolean", BCF_HT_FLAG},
      {"flag", BCF_HT_FLAG},
      {"Flag", BCF_HT_FLAG},
      {"string", BCF_HT_STR},
      {"String", BCF_HT_STR},
      {"char", BCF_HT_STR},
      {"Char", BCF_HT_STR}
      });

std::unordered_map<std::string, int> VidMapper::m_INFO_field_operation_name_to_enum =
  std::unordered_map<std::string, int>({
      {"sum", VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_SUM},
      {"mean", VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MEAN},
      {"median", VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MEDIAN},
      {"move_to_FORMAT", VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_MOVE_TO_FORMAT},
      {"element_wise_sum", VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_ELEMENT_WISE_SUM},
      {"concatenate", VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_CONCATENATE}
      });

#define VERIFY_OR_THROW(X) if(!(X)) throw VidMapperException(#X);

bool FileInfo::add_local_tiledb_row_idx_pair(const int64_t local, const int64_t global,
    int64_t& other_row_idx)
{
  auto iter = m_local_idx_to_tiledb_row_idx.find(local);
  //Already exists
  if(iter != m_local_idx_to_tiledb_row_idx.end())
  {
    other_row_idx = (*iter).second;
    //Given sample has been assigned two different row idxs - error
    if(other_row_idx != global)
      return false;
  }
  else
  {
    m_local_idx_to_tiledb_row_idx[local] = global;
    m_local_tiledb_row_idx_pairs.push_back(std::make_pair(local, global));
  }
  return true;
}

size_t FileInfo::get_num_orders() const
{
  size_t num_orders = 0u;
  switch(m_type)
  {
    case VidFileTypeEnum::VCF_FILE_TYPE:
    case VidFileTypeEnum::VCF_BUFFER_STREAM_TYPE:
    case VidFileTypeEnum::BCF_BUFFER_STREAM_TYPE:
    case VidFileTypeEnum::SORTED_CSV_FILE_TYPE:
    case VidFileTypeEnum::UNSORTED_CSV_FILE_TYPE:
      num_orders = 1u;
      break;
    default:
      throw VidMapperException(std::string("Unknown file type ")+std::to_string(m_type));
  }
  return num_orders;
}

void FieldLengthDescriptor::set_length_descriptor(const int length_dim_idx, const int length_descriptor)
{
  assert(static_cast<size_t>(length_dim_idx) < m_length_descriptor_vec.size());
  auto& curr_length_descriptor_component = m_length_descriptor_vec[length_dim_idx];
  curr_length_descriptor_component.m_length_descriptor = length_descriptor;
  m_is_fixed_length_field = m_is_fixed_length_field && (length_descriptor == BCF_VL_FIXED);
  m_is_length_genotype_dependent = m_is_length_genotype_dependent
    || KnownFieldInfo::is_length_descriptor_genotype_dependent(length_descriptor);
  m_is_length_allele_dependent = m_is_length_allele_dependent
    || KnownFieldInfo::is_length_descriptor_allele_dependent(length_descriptor);
  m_is_length_all_alleles_dependent = m_is_length_all_alleles_dependent
    || KnownFieldInfo::is_length_descriptor_all_alleles_dependent(length_descriptor);
  m_is_length_ploidy_dependent = m_is_length_ploidy_dependent
    || KnownFieldInfo::is_length_descriptor_ploidy_dependent(length_descriptor);
}

size_t FieldLengthDescriptor::get_num_elements(const unsigned num_ALT_alleles, const unsigned ploidy, const unsigned num_elements)
{
  assert(get_num_dimensions() == 1u);
  return KnownFieldInfo::get_num_elements_given_length_descriptor(get_length_descriptor(0u),
      num_ALT_alleles, ploidy, num_elements);
}

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
  m_buffer_stream_idx_to_global_file_idx.clear();
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
    if(field_info.m_is_vcf_FILTER_field)
      vcf_fields[BCF_HL_FLT].push_back(field_info.m_vcf_name);
    if(field_info.m_is_vcf_INFO_field)
      vcf_fields[BCF_HL_INFO].push_back(field_info.m_vcf_name);
    if(field_info.m_is_vcf_FORMAT_field)
      vcf_fields[BCF_HL_FMT].push_back(field_info.m_vcf_name);
  }
}

void VidMapper::build_tiledb_array_schema(VariantArraySchema*& array_schema, const std::string array_name,
    const bool row_based_partitioning, const RowRange& row_range, const bool compress_fields,
    const bool no_mandatory_VCF_fields)
  const
{
  auto dim_names = std::vector<std::string>({"samples", "position"});
  auto dim_domains = std::vector<std::pair<int64_t, int64_t>>({ {row_range.first, row_range.second}, {0, INT64_MAX}});
  std::vector<std::string> attribute_names;
  std::vector<std::type_index> types;
  std::vector<int> num_vals;
  //END
  attribute_names.push_back("END");
  types.push_back(std::type_index(typeid(int64_t)));
  num_vals.push_back(1);
  if(!no_mandatory_VCF_fields)
  {
    //REF
    attribute_names.push_back("REF");
    types.push_back(std::type_index(typeid(char)));
    num_vals.push_back(TILEDB_VAR_NUM);
    //ALT
    attribute_names.push_back("ALT");
    types.push_back(std::type_index(typeid(char)));
    num_vals.push_back(TILEDB_VAR_NUM);
    //ID field - optional
    auto iter = m_field_name_to_idx.find("ID");
    if(iter != m_field_name_to_idx.end())
    {
      attribute_names.push_back("ID");
      types.push_back(std::type_index(typeid(char)));
      num_vals.push_back(TILEDB_VAR_NUM);
    }
    //QUAL
    attribute_names.push_back("QUAL");
    types.push_back(std::type_index(typeid(float)));
    num_vals.push_back(1);
    //FILTER
    attribute_names.push_back("FILTER");
    types.push_back(std::type_index(typeid(int)));
    num_vals.push_back(TILEDB_VAR_NUM);
  }
  //INFO fields
  for(const auto& field_info : m_field_idx_to_info)
  {
    if(field_info.m_name == "END")      //skip END field
      continue;
    if(field_info.m_is_vcf_INFO_field)
    {
      attribute_names.push_back(field_info.m_name);
      types.push_back(field_info.m_type_index);
      num_vals.push_back(field_info.m_length_descriptor.is_fixed_length_field()
          ? field_info.m_length_descriptor.get_num_elements()
          : TILEDB_VAR_NUM);
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
      num_vals.push_back(field_info.m_length_descriptor.is_fixed_length_field()
          ? field_info.m_length_descriptor.get_num_elements()
          : TILEDB_VAR_NUM);
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
    //Either the file is not handled by the loader in this invocation or it's within the vector size
    VERIFY_OR_THROW(file_idx < 0 || static_cast<size_t>(file_idx) < m_file_idx_to_info.size());
    if(file_idx >= 0 && files_set.find(file_idx) == files_set.end())
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
    auto& curr_file_info = m_file_idx_to_info[global_file_idx];
    curr_file_info.m_local_file_idx = i;
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

std::string VidMapper::get_split_file_path(const std::string& original_filename, const std::string& results_directory,
    std::string& output_type, const int rank) const
{
  std::string return_value;
  auto iter = m_filename_to_idx.find(original_filename);
  //See if the callsets JSON had a path to output
  if(iter != m_filename_to_idx.end())
  {
    assert(static_cast<size_t>((*iter).second) < m_file_idx_to_info.size());
    const auto& file_info = m_file_idx_to_info[(*iter).second];
    if(file_info.m_single_split_file_path)
    {
      assert(file_info.m_split_files_paths.size() == 1u);
      return_value = file_info.m_split_files_paths[0u];
    }
    else
    {
      if(file_info.m_split_files_paths.size()) //split files were specified
      {
        if(static_cast<size_t>(rank) >= file_info.m_split_files_paths.size())
          throw VidMapperException(std::string("Split files entry for file ")+original_filename+" had "
                +std::to_string(file_info.m_split_files_paths.size())+" entries in the callset mapping JSON;"
                +" however, entry with index "+std::to_string(rank)+" is requested");
        return_value = file_info.m_split_files_paths[rank];
      }
    }
  }
  //No output filename specified - construct one in the same directory as the original file
  //Callsets JSON had no entry - create filepath
  if(return_value.empty())
  {
    //dirname and basename modify the original string - so pass copies
    auto copy_filepath = strdup(original_filename.c_str());
    std::string original_dirname = dirname(copy_filepath);
    //restore original
    memcpy(copy_filepath, original_filename.c_str(), original_filename.length());
    std::string original_basename = basename(copy_filepath);
    free(copy_filepath);
    return_value = ((results_directory.empty()) ? original_dirname : results_directory) + '/'
      + "partition_"+std::to_string(rank)+'_' + original_basename;
  }
  if(output_type.empty()) //deduce type by trying to open the original
  {
    auto fptr = bcf_open(original_filename.c_str(), "r");
    if(fptr)
    {
      switch(fptr->format.format)
      {
        case htsExactFormat::bcf:
          output_type = "b";
          break;
        case htsExactFormat::vcf:
          output_type = "z";
          break;
        default:
          break; //do nothing
      }
      bcf_close(fptr);
    }
  }
  return return_value;
}

std::vector<ContigIntervalTuple> VidMapper::get_contig_intervals_for_column_partition(
        const int64_t column_partition_begin, const int64_t column_partition_end, const bool is_zero_based) const
{
  auto to_add = is_zero_based ? 0 : 1;
  std::string contig_name = "";
  auto current_column = column_partition_begin;
  int64_t contig_position = -1ll;
  auto status = get_contig_location(current_column, contig_name, contig_position);
  std::vector<ContigIntervalTuple> contig_intervals;
  ContigInfo contig_info;
  if(status)
  {
    status = get_contig_info(contig_name, contig_info);
    assert(status);
    contig_intervals.emplace_back(contig_name, contig_position+to_add,
        contig_info.m_tiledb_column_offset+contig_info.m_length > column_partition_end //contig crosses the end of partition
        ? column_partition_end-contig_info.m_tiledb_column_offset+to_add
        : contig_info.m_length-1+to_add);
  }
  int64_t contig_offset = -1ll;
  auto next_contig_exists = get_next_contig_location(current_column, contig_name, contig_offset);
  while(next_contig_exists && contig_offset <= column_partition_end)
  {
    status = get_contig_info(contig_name, contig_info);
    assert(status);
    contig_intervals.emplace_back(contig_name, to_add,
        contig_info.m_tiledb_column_offset+contig_info.m_length > column_partition_end //contig crosses the end of partition
        ? column_partition_end-contig_info.m_tiledb_column_offset+to_add
        : contig_info.m_length-1+to_add);
    current_column = contig_offset;
    next_contig_exists = get_next_contig_location(current_column, contig_name, contig_offset);
  }
  return contig_intervals;
}

std::vector<ContigIntervalTuple> VidMapper::get_contig_intervals_for_column_partition(
        const std::string& loader_filename,
        const int rank, const bool is_zero_based)
{
  JSONLoaderConfig loader_config;
  loader_config.read_from_file(loader_filename);
  //don't need a callset file
  FileBasedVidMapper vid_mapper(loader_config.get_vid_mapping_filename(), "", 0, INT64_MAX-1, false);
  if(loader_config.is_partitioned_by_row())
    return vid_mapper.get_contig_intervals_for_column_partition(0, INT64_MAX-1, is_zero_based);
  else
  {
    auto column_partition = loader_config.get_column_partition(rank);
    return vid_mapper.get_contig_intervals_for_column_partition(column_partition.first, column_partition.second, is_zero_based);
  }
}

//FileBasedVidMapper code
#ifdef VERIFY_OR_THROW
#undef VERIFY_OR_THROW
#endif

#define VERIFY_OR_THROW(X) if(!(X)) throw FileBasedVidMapperException(#X);

void FileBasedVidMapper::common_constructor_initialization(const std::string& filename,
    const std::vector<BufferStreamInfo>& buffer_stream_info_vec,
    const std::string& callset_mapping_file,
    const std::string& buffer_stream_callset_mapping_json_string,
    const int64_t lb_callset_row_idx, const int64_t ub_callset_row_idx, const bool is_callset_mapping_required)
{
  m_lb_callset_row_idx = 0;
  m_ub_callset_row_idx = INT64_MAX-1;
  VERIFY_OR_THROW(filename.length() && "Vid mapping file unspecified");
  rapidjson::Document json_doc;
  std::ifstream ifs(filename.c_str());
  if(!ifs.is_open())
    throw FileBasedVidMapperException((std::string("Could not open vid mapping file \"")+filename+"\"").c_str());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  json_doc.Parse(str.c_str());
  if(json_doc.HasParseError())
    throw FileBasedVidMapperException(std::string("Syntax error in JSON file ")+filename);
  m_lb_callset_row_idx = lb_callset_row_idx;
  m_ub_callset_row_idx = ub_callset_row_idx;
  //Callset info parsing
  std::string real_callset_mapping_file;
  if(!callset_mapping_file.empty())
    real_callset_mapping_file = callset_mapping_file;
  else
  {
    if(json_doc.HasMember("callset_mapping_file"))
    {
      VERIFY_OR_THROW(json_doc["callset_mapping_file"].IsString() && "The \"callset_mapping_file\" field must be a string and specify a path to the callsets file");
      real_callset_mapping_file = json_doc["callset_mapping_file"].GetString();
    }
    else
      if(buffer_stream_callset_mapping_json_string.empty())
        if(is_callset_mapping_required)
          throw FileBasedVidMapperException(std::string("Must specify callset mapping either through \"callset_mapping_file\" or through a JSON string for buffer streams"));
  }
  parse_callsets_json(real_callset_mapping_file, buffer_stream_info_vec, true);
  parse_callsets_json(buffer_stream_callset_mapping_json_string, buffer_stream_info_vec, false);
  //Contig info parsing
  VERIFY_OR_THROW(json_doc.HasMember("contigs"));
  {
    const rapidjson::Value& contigs_container = json_doc["contigs"];
    //contigs is a dictionary of name:info key-value pairs or array of dict { "name": <>, ... }
    VERIFY_OR_THROW(contigs_container.IsObject() || contigs_container.IsArray());
    auto is_array = contigs_container.IsArray();
    auto num_contigs = is_array ? contigs_container.Size() : contigs_container.MemberCount();
    m_contig_idx_to_info.resize(num_contigs);
    m_contig_begin_2_idx.resize(num_contigs);
    m_contig_end_2_idx.resize(num_contigs);
    std::string contig_name;
    std::string filename;
    auto duplicate_contigs_exist = false;
    //index within the JSON array
    rapidjson::SizeType json_contig_idx = 0u;
    auto dict_iter = is_array ? rapidjson::Value::ConstMemberIterator()
      : contigs_container.MemberBegin();
    auto dict_end_position = is_array ? dict_iter : contigs_container.MemberEnd();
    auto next_contig_exists = is_array ? (json_contig_idx < num_contigs) : (dict_iter != dict_end_position);
    while(next_contig_exists)
    {
      const auto& contig_info_dict = is_array ? contigs_container[json_contig_idx] : (*dict_iter).value;
      VERIFY_OR_THROW(contig_info_dict.IsObject());
      if(is_array)
      {
        auto num_found = 0u;
        for(const auto name_field : { "name", "contig_name", "chromosome_name" })
        {
          if(contig_info_dict.HasMember(name_field))
          {
            contig_name = contig_info_dict[name_field].GetString();
            ++num_found;
          }
        }
        if(num_found == 0u)
          throw VidMapperException(std::string("Contig info dict with index ")+std::to_string(json_contig_idx)
              +" does not have any one of the keys \"name\", \"contig_name\" or \"chromosome_name\"");
        if(num_found > 1u)
          throw VidMapperException(std::string("Contig info dict with index ")+std::to_string(json_contig_idx)
              +" has two or more of the keys \"name\", \"contig_name\" or \"chromosome_name\" - at most one is allowed");
      }
      else
        contig_name = (*dict_iter).name.GetString();
      if(m_contig_name_to_idx.find(contig_name) != m_contig_name_to_idx.end())
      {
        std::cerr << (std::string("Duplicate contig/chromosome name ")+contig_name+" found in vid file "+filename) << "\n";
        duplicate_contigs_exist = true;
      }
      else
      {
        VERIFY_OR_THROW(contig_info_dict.HasMember("tiledb_column_offset") && contig_info_dict["tiledb_column_offset"].IsInt64());
        auto tiledb_column_offset = contig_info_dict["tiledb_column_offset"].GetInt64();
        VERIFY_OR_THROW(tiledb_column_offset >= 0ll);
        VERIFY_OR_THROW(contig_info_dict.HasMember("length") && contig_info_dict["length"].IsInt64());
        auto length = contig_info_dict["length"].GetInt64();
        VERIFY_OR_THROW(length >= 0ll);
        VERIFY_OR_THROW(static_cast<size_t>(json_contig_idx) < static_cast<size_t>(num_contigs));
        m_contig_name_to_idx[contig_name] = json_contig_idx;
        m_contig_idx_to_info[json_contig_idx].set_info(json_contig_idx, contig_name, length, tiledb_column_offset);
        m_contig_begin_2_idx[json_contig_idx].first = tiledb_column_offset;
        m_contig_begin_2_idx[json_contig_idx].second = json_contig_idx;
        m_contig_end_2_idx[json_contig_idx].first = tiledb_column_offset + length - 1; //inclusive
        m_contig_end_2_idx[json_contig_idx].second = json_contig_idx;
      }
      ++json_contig_idx;
      if(!is_array)
        ++dict_iter;
      next_contig_exists = is_array ? (json_contig_idx < num_contigs) : (dict_iter != dict_end_position);
    }
    if(duplicate_contigs_exist)
      throw FileBasedVidMapperException(std::string("Duplicate contigs exist in vid file ")+filename);
    std::sort(m_contig_begin_2_idx.begin(), m_contig_begin_2_idx.end(), contig_offset_idx_pair_cmp);
    std::sort(m_contig_end_2_idx.begin(), m_contig_end_2_idx.end(), contig_offset_idx_pair_cmp);
    //Check that there are no spurious overlaps
    auto last_contig_idx = -1;
    auto last_contig_end_column = -1ll;
    auto overlapping_contigs_exist = false;
    for(const auto& offset_idx_pair : m_contig_begin_2_idx)
    {
      auto contig_idx = offset_idx_pair.second;
      const auto& contig_info = m_contig_idx_to_info[contig_idx];
      if(last_contig_idx >= 0)
      {
        const auto& last_contig_info = m_contig_idx_to_info[last_contig_idx];
        if(contig_info.m_tiledb_column_offset <= last_contig_end_column)
        {
          std::cerr << (std::string("Contig/chromosome ")+contig_info.m_name+" begins at TileDB column "
              +std::to_string(contig_info.m_tiledb_column_offset)+" and intersects with contig/chromosome "+last_contig_info.m_name
              +" that spans columns [ "+std::to_string(last_contig_info.m_tiledb_column_offset)+", "
              +std::to_string(last_contig_info.m_tiledb_column_offset+last_contig_info.m_length-1)+" ]") << "\n";
          overlapping_contigs_exist = true;
        }
      }
      last_contig_idx = contig_idx;
      last_contig_end_column = contig_info.m_tiledb_column_offset + contig_info.m_length - 1;
    }
    if(overlapping_contigs_exist)
      throw FileBasedVidMapperException(std::string("Overlapping contigs exist in vid file ")+filename);
  }
  //Field info parsing
  VERIFY_OR_THROW(json_doc.HasMember("fields"));
  {
    const rapidjson::Value& fields_container = json_doc["fields"];
    //fields is a dictionary: name: dict or array[dict]
    VERIFY_OR_THROW(fields_container.IsObject() || fields_container.IsArray());
    auto is_array = fields_container.IsArray();
    auto num_fields = is_array ? fields_container.Size() : fields_container.MemberCount();
    m_field_idx_to_info.resize(num_fields);
    std::string field_name;
    auto duplicate_fields_exist = false;
    //index within m_field_idx_to_info
    auto field_idx = 0;
    //index within the JSON array
    rapidjson::SizeType json_field_idx = 0u;
    auto dict_iter = is_array ? rapidjson::Value::ConstMemberIterator()
      : fields_container.MemberBegin();
    auto dict_end_position = is_array ? dict_iter : fields_container.MemberEnd();
    auto next_field_exists = is_array ? (json_field_idx < fields_container.Size()) : (dict_iter != dict_end_position);
    while(next_field_exists)
    {
      const auto& field_info_dict = is_array ? fields_container[json_field_idx] : (*dict_iter).value;
      VERIFY_OR_THROW(field_info_dict.IsObject());
      if(is_array)
      {
        if(!field_info_dict.HasMember("name") && !field_info_dict.HasMember("field_name"))
          throw VidMapperException(std::string("Field information dictionary with index ")+std::to_string(json_field_idx)
                +" does not have a \"field_name\" or \"name\" key");
        if(field_info_dict.HasMember("name") && field_info_dict.HasMember("field_name"))
          throw VidMapperException(std::string("Field information dictionary with index ")+std::to_string(json_field_idx)
                +" has both \"field_name\" or \"name\" keys - only one should be specified");
        field_name = field_info_dict.HasMember("name") ? field_info_dict["name"].GetString()
          : field_info_dict["field_name"].GetString();
      }
      else
        field_name = (*dict_iter).name.GetString();
      if(m_field_name_to_idx.find(field_name) != m_field_name_to_idx.end())
      {
        std::cerr << (std::string("Duplicate field name ")+field_name+" found in vid file "+filename) << "\n";
        duplicate_fields_exist = true;
        //advance loop
      }
      else
      {
        //Known fields
        auto known_field_enum = 0u;
        auto is_known_field = KnownFieldInfo::get_known_field_enum_for_name(field_name, known_field_enum);
        //Map
        m_field_name_to_idx[field_name] = field_idx;
        m_field_idx_to_info[field_idx].set_info(field_name, field_idx);
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
              else
                if(class_name == "FILTER")
                  m_field_idx_to_info[field_idx].m_is_vcf_FILTER_field = true;
          }
        }
        if(field_info_dict.HasMember("length"))
          parse_length_descriptor(field_name.c_str(),
              field_info_dict["length"], m_field_idx_to_info[field_idx].m_length_descriptor, 0u);
        else
        {
          if(is_known_field)
          {
            auto length_descriptor_code = KnownFieldInfo::get_length_descriptor_for_known_field_enum(known_field_enum);
            m_field_idx_to_info[field_idx].m_length_descriptor.set_length_descriptor(0u, length_descriptor_code);
            if(length_descriptor_code == BCF_VL_FIXED)
              m_field_idx_to_info[field_idx].m_length_descriptor.set_num_elements(0u,
                  KnownFieldInfo::get_num_elements_for_known_field_enum(known_field_enum, 0u, 0u));  //don't care about ploidy
          }
        }
        if(field_info_dict.HasMember("VCF_field_combine_operation"))
        {
          VERIFY_OR_THROW(field_info_dict["VCF_field_combine_operation"].IsString());
          auto iter  = VidMapper::m_INFO_field_operation_name_to_enum.find(field_info_dict["VCF_field_combine_operation"].GetString());
          if(iter == VidMapper::m_INFO_field_operation_name_to_enum.end())
            throw VidMapperException(std::string("Unknown VCF field combine operation ")+field_info_dict["VCF_field_combine_operation"].GetString()
                +" specified for field "+field_name);
          m_field_idx_to_info[field_idx].m_VCF_field_combine_operation = (*iter).second;
          //Concatenate can only be used for VAR length fields
          if(m_field_idx_to_info[field_idx].m_VCF_field_combine_operation == VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_CONCATENATE
              && m_field_idx_to_info[field_idx].m_length_descriptor.get_length_descriptor(0u) != BCF_VL_VAR)
            throw VidMapperException(std::string("VCF field combined operation 'concatenate' can only be used with fields whose length descriptors are 'VAR'; ")
                +" field "+field_name+" does not have 'VAR' as its length descriptor");
        }
        else
        {
          if(is_known_field)
            m_field_idx_to_info[field_idx].m_VCF_field_combine_operation = KnownFieldInfo::get_VCF_field_combine_operation_for_known_field_enum(known_field_enum);
        }
        //Both INFO and FORMAT, throw another entry <field>_FORMAT
        if(m_field_idx_to_info[field_idx].m_is_vcf_INFO_field && m_field_idx_to_info[field_idx].m_is_vcf_FORMAT_field)
        {
          auto new_field_idx = field_idx+1u;
          m_field_idx_to_info.resize(m_field_idx_to_info.size()+1u);
          //Copy field information
          m_field_idx_to_info[new_field_idx] = m_field_idx_to_info[field_idx];
          auto& new_field_info =  m_field_idx_to_info[new_field_idx];
          //Update name and index - keep the same VCF name
          new_field_info.m_name = field_name+"_FORMAT";
          new_field_info.m_is_vcf_INFO_field = false;
          new_field_info.m_field_idx = new_field_idx;
          new_field_info.m_VCF_field_combine_operation = VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_UNKNOWN_OPERATION;
          //Update map
          m_field_name_to_idx[new_field_info.m_name] = new_field_idx;
          //Set FORMAT to false for original field
          m_field_idx_to_info[field_idx].m_is_vcf_FORMAT_field = false;
          ++field_idx;
        }
      }
      ++field_idx;
      ++json_field_idx;
      if(!is_array)
        ++dict_iter;
      next_field_exists = is_array ? (json_field_idx < fields_container.Size()) : (dict_iter != dict_end_position);
    }
    //Force add END as a field
    auto iter = m_field_name_to_idx.find("END");
    if(iter == m_field_name_to_idx.end())
    {
      auto end_idx = m_field_idx_to_info.size();
      m_field_idx_to_info.emplace_back();
      m_field_name_to_idx["END"] = end_idx;
      auto& field_info = m_field_idx_to_info[end_idx];
      field_info.set_info("END", end_idx);
      field_info.m_is_vcf_INFO_field = true;
      field_info.m_type_index = std::move(std::type_index(typeid(int)));
      field_info.m_bcf_ht_type = BCF_HT_INT;
    }
    if(duplicate_fields_exist)
      throw FileBasedVidMapperException(std::string("Duplicate fields exist in vid file ")+filename);
  }
  m_is_initialized = true;
}

void FileBasedVidMapper::parse_length_descriptor(const char* field_name,
    const rapidjson::Value& length_json_value,
    FieldLengthDescriptor& length_descriptor, const size_t length_dim_idx)
{
  if(length_json_value.IsInt64())
    length_descriptor.set_num_elements(length_dim_idx, length_json_value.GetInt64());
  else
    if(length_json_value.IsString())
      parse_string_length_descriptor(field_name,
          length_json_value.GetString(),
          length_json_value.GetStringLength(),
          length_descriptor, length_dim_idx);
    else
    {
      //Protobuf produces a JSON which looks like this:
      //"length" : [ { "variable_length_descriptor": "R" },  { "fixed_length" : 2 } ]
      if(length_json_value.IsObject())
      {
        if(length_json_value.HasMember("variable_length_descriptor"))
        {
          auto& str_value = length_json_value["variable_length_descriptor"];
          VERIFY_OR_THROW(str_value.IsString());
          parse_string_length_descriptor(field_name,
              str_value.GetString(),
              str_value.GetStringLength(),
              length_descriptor, length_dim_idx);
        }
        else
        {
          VERIFY_OR_THROW(length_json_value.HasMember("fixed_length"));
          auto& int_value = length_json_value["fixed_length"];
          VERIFY_OR_THROW(int_value.IsInt64());
          length_descriptor.set_num_elements(length_dim_idx, int_value.GetInt64());
        }
      }
      else
      {
        //MultiD array of form [ "R", "var", { "variable_length_descriptor": "A" } ] ..
        VERIFY_OR_THROW(length_json_value.IsArray());
        length_descriptor.resize(length_json_value.Size());
        for(rapidjson::SizeType idx=0u;idx<length_json_value.Size();++idx)
          parse_length_descriptor(field_name, length_json_value[idx], length_descriptor, idx);
      }
    }
}

void VidMapper::parse_string_length_descriptor(
    const char* field_name,
    const char* length_value_str,
    const size_t length_value_str_length,
    FieldLengthDescriptor& length_descriptor, const size_t length_dim_idx)
{
  auto length_value_upper_case_str = std::move(std::string(length_value_str));
  for(auto i=0u;i<length_value_upper_case_str.length();++i)
    length_value_upper_case_str[i] = toupper(length_value_upper_case_str[i]);
  auto iter = VidMapper::m_length_descriptor_string_to_int.find(length_value_upper_case_str);
  if(iter == VidMapper::m_length_descriptor_string_to_int.end())
  {
    //JSON produced by Protobuf specifies fixed length field lengths as strings - e.g. "1"
    char* endptr = 0;
    auto length_value_int = strtoull(length_value_str, &endptr, 0);
    auto num_chars_traversed = endptr-length_value_str;
    if(length_value_str_length > 0u
        && static_cast<size_t>(num_chars_traversed) == length_value_str_length) //whole string is an integer
      length_descriptor.set_num_elements(length_dim_idx, length_value_int);
    else
    {
      std::cerr << "WARNING: unknown length descriptor " << length_value_str
        << " for field " << field_name  << " ; setting to 'VAR'\n";
      length_descriptor.set_length_descriptor(length_dim_idx, BCF_VL_VAR);
    }
  }
  else
    length_descriptor.set_length_descriptor(length_dim_idx, (*iter).second);
}

void FileBasedVidMapper::parse_callsets_json(const std::string& json, const std::vector<BufferStreamInfo>& buffer_stream_info_vec,
    const bool is_file)
{
  if(json.empty())
    return;
  std::string filename = is_file ? json : "buffer_stream_callset_mapping_json_string";
  rapidjson::Document json_doc;
  if(is_file)
  {
    std::ifstream ifs(filename.c_str());
    if(!ifs.is_open())
      throw FileBasedVidMapperException((std::string("Could not open callsets file \"")+filename+"\"").c_str());
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    json_doc.Parse(str.c_str());
    if(json_doc.HasParseError())
      throw FileBasedVidMapperException(std::string("Syntax error in JSON file ")+filename);
  }
  else
  {
    json_doc.Parse(json.c_str());
    if(json_doc.HasParseError())
      throw FileBasedVidMapperException(std::string("Syntax error in JSON ")+filename);
  }
  //Callset info parsing
  VERIFY_OR_THROW(json_doc.HasMember("callsets"));
  {
    const rapidjson::Value& callsets_container = json_doc["callsets"];
    //callsets is a dictionary of name:info key-value pairs or array of info dictionaries
    VERIFY_OR_THROW(callsets_container.IsObject() || callsets_container.IsArray());
    auto is_array = callsets_container.IsArray();
    auto num_callsets_in_container = is_array ? callsets_container.Size() : callsets_container.MemberCount();
    auto num_callsets = (m_ub_callset_row_idx == INT64_MAX) ? num_callsets_in_container :
      std::min<int64_t>(num_callsets_in_container, m_ub_callset_row_idx+1);
    m_row_idx_to_info.resize(num_callsets);
    m_max_callset_row_idx = -1;
    std::string callset_name;
    auto dict_iter = is_array ? rapidjson::Value::ConstMemberIterator()
      : callsets_container.MemberBegin();
    auto dict_end_position = is_array ? dict_iter : callsets_container.MemberEnd();
    rapidjson::SizeType json_callset_idx = 0ull;
    auto next_callset_exists = is_array ? (json_callset_idx < callsets_container.Size()) : (dict_iter != dict_end_position);
    while(next_callset_exists)
    {
      const auto& callset_info_dict = is_array ? callsets_container[json_callset_idx]
        : (*dict_iter).value;
      VERIFY_OR_THROW(callset_info_dict.IsObject());    //must be dict
      if(is_array)
      {
        auto num_found = 0u;
        for(const auto name_field : { "name", "callset_name", "sample_name" })
        {
          if(callset_info_dict.HasMember(name_field))
          {
            callset_name = callset_info_dict[name_field].GetString();
            ++num_found;
          }
        }
        if(num_found == 0u)
          throw VidMapperException(std::string("Callset info dict with index ")+std::to_string(json_callset_idx)
              +" does not have any one of the keys \"name\", \"callset_name\" or \"sample_name\"");
        if(num_found > 1u)
          throw VidMapperException(std::string("Callset info dict with index ")+std::to_string(json_callset_idx)
              +" has two or more of the keys \"name\", \"callset_name\" or \"sample_name\" - at most one is allowed");
      }
      else
        callset_name = (*dict_iter).name.GetString();
      VERIFY_OR_THROW(callset_info_dict.HasMember("row_idx"));
      int64_t row_idx = callset_info_dict["row_idx"].GetInt64();
      //already exists in map
      if(m_callset_name_to_row_idx.find(callset_name) != m_callset_name_to_row_idx.end())
        //different row idx
        if(m_callset_name_to_row_idx[callset_name] != row_idx)
          throw FileBasedVidMapperException(std::string("Duplicate with conflicting row index for sample/callset name ")+callset_name
              +" found in callsets mapping "+std::to_string(m_callset_name_to_row_idx[callset_name])+", "+std::to_string(row_idx));
      if(row_idx <= m_ub_callset_row_idx)
      {
        m_max_callset_row_idx = std::max(m_max_callset_row_idx, row_idx);
        //Resize vector
        if(static_cast<size_t>(row_idx) >= m_row_idx_to_info.size())
          m_row_idx_to_info.resize(row_idx+1);
        auto file_idx = -1ll;
        //idx in file
        auto idx_in_file = 0ll;
        if(row_idx >= m_lb_callset_row_idx && (callset_info_dict.HasMember("filename") || callset_info_dict.HasMember("stream_name")))
        {
          VERIFY_OR_THROW((!callset_info_dict.HasMember("filename") || !callset_info_dict.HasMember("stream_name"))
              && (std::string("Cannot have both \"filename\" and \"stream_name\" as the data source for sample/CallSet ")+callset_name).c_str());
          std::string filename = callset_info_dict.HasMember("filename")
            ? std::move(callset_info_dict["filename"].GetString())
            : std::move(callset_info_dict["stream_name"].GetString());
          file_idx = get_or_append_global_file_idx(filename);
          if(callset_info_dict.HasMember("idx_in_file"))
            idx_in_file = callset_info_dict["idx_in_file"].GetInt64();
          //Check for conflicting file/stream info if initialized previously
          const auto& curr_row_info = m_row_idx_to_info[row_idx];
          if(curr_row_info.m_is_initialized)
          {
            if(curr_row_info.m_file_idx >= 0 && file_idx >= 0
                && curr_row_info.m_file_idx != file_idx)
              throw FileBasedVidMapperException(std::string("Conflicting file/stream names specified for sample/callset ")+callset_name
                  +" "+m_file_idx_to_info[curr_row_info.m_file_idx].m_name+", "+filename);
            if(curr_row_info.m_idx_in_file != idx_in_file)
              throw FileBasedVidMapperException(std::string("Conflicting values of \"idx_in_file\" specified for sample/callset ")+callset_name
                  +" "+std::to_string(curr_row_info.m_idx_in_file)+", "+std::to_string(idx_in_file));
          }
          else
          {
            assert(file_idx < static_cast<int64_t>(m_file_idx_to_info.size()));
            int64_t other_row_idx = 0ll;
            auto added_successfully = m_file_idx_to_info[file_idx].add_local_tiledb_row_idx_pair(idx_in_file,
                row_idx, other_row_idx);
            if(!added_successfully)
              throw FileBasedVidMapperException(std::string("Attempting to import a sample from file/stream ")+filename
                  +" multiple times under aliases '"+m_row_idx_to_info[other_row_idx].m_name
                  +"' and '"+callset_name+"' with row indexes "+std::to_string(other_row_idx)
                  +" and "+std::to_string(row_idx)+" respectively");
          }
        }
        m_callset_name_to_row_idx[callset_name] = row_idx;
        VERIFY_OR_THROW(static_cast<size_t>(row_idx) < m_row_idx_to_info.size());
        if(m_row_idx_to_info[row_idx].m_is_initialized && m_row_idx_to_info[row_idx].m_name != callset_name)
          throw FileBasedVidMapperException(std::string("Sample/callset ")+callset_name+" has the same row idx as "
              +m_row_idx_to_info[row_idx].m_name);
        m_row_idx_to_info[row_idx].set_info(row_idx, callset_name, file_idx, idx_in_file);
      }
      ++json_callset_idx;
      if(!is_array)
        ++dict_iter;
      next_callset_exists = is_array ? (json_callset_idx < callsets_container.Size()) : (dict_iter != dict_end_position);
    }
    auto missing_row_idxs_exist = false;
    for(auto row_idx=0ll;static_cast<size_t>(row_idx)<m_row_idx_to_info.size() && row_idx<=m_ub_callset_row_idx;++row_idx)
      if(row_idx >= m_lb_callset_row_idx && row_idx <= m_ub_callset_row_idx && !(m_row_idx_to_info[row_idx].m_is_initialized))
      {
        std::cerr << "Sample/callset information missing for row " << row_idx << "\n";
        missing_row_idxs_exist = true;
      }
    if(missing_row_idxs_exist)
      throw FileBasedVidMapperException(std::string("Row indexes with missing sample/callset information found in callsets mapping: ")+filename);
  }
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
  //File/stream types
  for(const auto& entry : std::unordered_map<std::string, VidFileTypeEnum>({
        {"sorted_csv_files", VidFileTypeEnum::SORTED_CSV_FILE_TYPE },
        {"unsorted_csv_files", VidFileTypeEnum::UNSORTED_CSV_FILE_TYPE },
        {"vcf_buffer_streams", VidFileTypeEnum::VCF_BUFFER_STREAM_TYPE },
        {"bcf_buffer_streams", VidFileTypeEnum::BCF_BUFFER_STREAM_TYPE }
        }))
  {
    if(json_doc.HasMember(entry.first.c_str()))
    {
      const auto& file_array = json_doc[entry.first.c_str()];
      //is an array
      VERIFY_OR_THROW(file_array.IsArray());
      for(rapidjson::SizeType i=0;i<file_array.Size();++i)
      {
        int64_t global_file_idx;
        auto found = get_global_file_idx(file_array[i].GetString(), global_file_idx);
        if(found)
          m_file_idx_to_info[global_file_idx].m_type = entry.second;
      }
    }
  }
  //If splitting VCFs, might have split file paths
  if(json_doc.HasMember("split_files_paths"))
  {
    //json_doc["split_files_paths"] must be dictionary of the form:
    //<original_file>: [ "split_paths" ] or
    //<original_file>: "split_path"
    VERIFY_OR_THROW(json_doc["split_files_paths"].IsObject());
    const auto& split_files_dict = json_doc["split_files_paths"];
    std::unordered_set<std::string> check_duplicates_set;
    for(auto b=split_files_dict.MemberBegin(), e=split_files_dict.MemberEnd();b!=e;++b)
    {
      const auto& curr_obj = *b;
      const auto& original_filename = curr_obj.name.GetString();
      if(check_duplicates_set.find(original_filename) != check_duplicates_set.end())
        throw FileBasedVidMapperException(std::string("File ")+original_filename+" listed multiple times inside \"split_files_paths\" in the callsets mapping file");
      check_duplicates_set.insert(original_filename);
      //Insert if missing
      auto file_idx = get_or_append_global_file_idx(original_filename);
      assert(static_cast<size_t>(file_idx) < m_file_idx_to_info.size());
      auto& file_info = m_file_idx_to_info[file_idx];
      const auto& split_files_values = curr_obj.value;
      VERIFY_OR_THROW(split_files_values.IsArray() || split_files_values.IsString());
      if(split_files_values.IsArray())
      {
        for(rapidjson::SizeType i=0u;i<split_files_values.Size();++i)
        {
          VERIFY_OR_THROW(split_files_values[i].IsString());
          file_info.m_split_files_paths.push_back(split_files_values[i].GetString());
        }
      }
      else
      {
        file_info.m_single_split_file_path = true;
        file_info.m_split_files_paths.push_back(split_files_values.GetString());
      }
    }
  }
  //For buffer streams
  m_buffer_stream_idx_to_global_file_idx.resize(buffer_stream_info_vec.size(), -1);
  auto max_buffer_stream_idx_with_global_file_idx = -1ll;
  for(auto i=0ull;i<buffer_stream_info_vec.size();++i)
  {
    const auto& info = buffer_stream_info_vec[i];
    int64_t global_file_idx;
    auto found = get_global_file_idx(info.m_name, global_file_idx);
    if(found)
    {
      auto& curr_file_info = m_file_idx_to_info[global_file_idx];
      curr_file_info.m_type = info.m_type;
      curr_file_info.m_buffer_stream_idx = i;
      curr_file_info.m_buffer_capacity = info.m_buffer_capacity;
      curr_file_info.m_initialization_buffer = info.m_initialization_buffer;
      curr_file_info.m_initialization_buffer_num_valid_bytes = info.m_initialization_buffer_num_valid_bytes;
      m_buffer_stream_idx_to_global_file_idx[i] = global_file_idx;
      max_buffer_stream_idx_with_global_file_idx = i;
    }
  }
  m_buffer_stream_idx_to_global_file_idx.resize(max_buffer_stream_idx_with_global_file_idx+1);
  m_is_callset_mapping_initialized = true;
}

void FileBasedVidMapper::write_partition_callsets_json_file(const std::string& original_callsets_filename, const std::string& results_directory,
    const int rank) const
{
  rapidjson::Document json_doc(rapidjson::kObjectType);
  json_doc.SetObject();
  auto& allocator = json_doc.GetAllocator();
  rapidjson::Value json_string_value(rapidjson::kStringType);
  //Callsets dictionary
  std::string output_type;
  rapidjson::Value callsets_dict(rapidjson::kObjectType);
  for(auto i=0ull;i<m_row_idx_to_info.size();++i)
  {
    const auto& curr_callset_info = m_row_idx_to_info[i];
    rapidjson::Value curr_callset_dict(rapidjson::kObjectType);
    curr_callset_dict.AddMember("row_idx", curr_callset_info.m_row_idx, allocator);
    curr_callset_dict.AddMember("idx_in_file", curr_callset_info.m_idx_in_file, allocator);
    if(curr_callset_info.m_file_idx >= 0)
    {
      assert(static_cast<size_t>(curr_callset_info.m_file_idx) < m_file_idx_to_info.size());
      auto& original_filename = m_file_idx_to_info[curr_callset_info.m_file_idx].m_name;
      output_type.clear();
      auto output_filename = get_split_file_path(original_filename, results_directory, output_type, rank);
      json_string_value.SetString(output_filename.c_str(), output_filename.length(), allocator);
      curr_callset_dict.AddMember("filename", json_string_value, allocator);
    }
    callsets_dict.AddMember(rapidjson::StringRef(curr_callset_info.m_name.c_str()), curr_callset_dict, allocator);
  }
  json_doc.AddMember("callsets", callsets_dict, allocator);
  //File/stream types
  for(const auto& entry : std::unordered_map<std::string, VidFileTypeEnum>({
        {"sorted_csv_files", VidFileTypeEnum::SORTED_CSV_FILE_TYPE },
        {"unsorted_csv_files", VidFileTypeEnum::UNSORTED_CSV_FILE_TYPE },
        {"vcf_buffer_streams", VidFileTypeEnum::VCF_BUFFER_STREAM_TYPE },
        {"bcf_buffer_streams", VidFileTypeEnum::BCF_BUFFER_STREAM_TYPE }
        }))
  {
    rapidjson::Value file_type_list(rapidjson::kArrayType);
    for(auto i=0ull;i<m_file_idx_to_info.size();++i)
    {
      auto& curr_file_info = m_file_idx_to_info[i];
      if(curr_file_info.m_type == entry.second)
      {
        output_type.clear();
        auto output_filename = get_split_file_path(curr_file_info.m_name, results_directory, output_type, rank);
        json_string_value.SetString(output_filename.c_str(), output_filename.length(), allocator);
        file_type_list.PushBack(json_string_value, allocator);
      }
    }
    if(file_type_list.Size() > 0u)
    {
      json_string_value.SetString(entry.first.c_str(), entry.first.length(), allocator);
      json_doc.AddMember(json_string_value, file_type_list, allocator);
    }
  }
  output_type.clear();
  auto output_filename = get_split_file_path(original_callsets_filename, results_directory, output_type, rank);
  auto* fptr = fopen(output_filename.c_str(), "w");
  if(fptr == 0)
    throw FileBasedVidMapperException(std::string("Could not write to partitioned callsets JSON file ")+output_filename);
  char write_buffer[65536];
  rapidjson::FileWriteStream writer_stream(fptr, write_buffer, sizeof(write_buffer));
  rapidjson::PrettyWriter<rapidjson::FileWriteStream> writer(writer_stream);
  json_doc.Accept(writer);
  fclose(fptr);
}

void FileBasedVidMapper::write_partition_loader_json_file(const std::string& original_loader_filename,
    const std::string& original_callsets_filename,
    const std::string& results_directory, const int num_callset_mapping_files, const int rank) const
{
  //Parse original loader json
  rapidjson::Document json_doc;
  std::ifstream ifs(original_loader_filename.c_str());
  VERIFY_OR_THROW(ifs.is_open());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  json_doc.Parse(str.c_str());
  if(json_doc.HasParseError())
    throw RunConfigException(std::string("Syntax error in JSON file ")+original_loader_filename);
  auto& allocator = json_doc.GetAllocator();
  //Partitioned callsets JSON
  std::string output_type;
  if(json_doc.HasMember("callset_mapping_file"))
    json_doc.RemoveMember("callset_mapping_file");
  rapidjson::Value json_string_value(rapidjson::kStringType);
  rapidjson::Value partitioned_callset_mapping_file_array(rapidjson::kArrayType);
  for(auto i=0;i<num_callset_mapping_files;++i)
  {
    auto output_filename = get_split_file_path(original_callsets_filename, results_directory, output_type,
        (num_callset_mapping_files > 1) ? i : rank);
    json_string_value.SetString(output_filename.c_str(), output_filename.length(), allocator);
    partitioned_callset_mapping_file_array.PushBack(json_string_value, allocator);
  }
  json_doc.AddMember("callset_mapping_file", (num_callset_mapping_files > 1) ? partitioned_callset_mapping_file_array : json_string_value,
      allocator);
  output_type.clear();
  auto output_filename = get_split_file_path(original_loader_filename, results_directory, output_type, rank);
  auto* fptr = fopen(output_filename.c_str(), "w");
  if(fptr == 0)
    throw FileBasedVidMapperException(std::string("Could not write to partitioned loader JSON file ")+output_filename);
  char write_buffer[65536];
  rapidjson::FileWriteStream writer_stream(fptr, write_buffer, sizeof(write_buffer));
  rapidjson::PrettyWriter<rapidjson::FileWriteStream> writer(writer_stream);
  json_doc.Accept(writer);
  fclose(fptr);
}
