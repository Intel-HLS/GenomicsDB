/**
 * The MIT License (MIT)
 * Copyright (c) 2018 Intel Corporation
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
#include "genomicsdb_config_base.h"
#include "tiledb_utils.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw GenomicsDBConfigException(#X);
    
std::unordered_map<std::string, bool> GenomicsDBConfigBase::m_vcf_output_format_to_is_bcf_flag =
std::unordered_map<std::string, bool>{ {"b", true}, {"bu",true}, {"z",false}, {"",false} };

ColumnRange GenomicsDBConfigBase::verify_contig_position_and_get_tiledb_column_interval(const ContigInfo& contig_info,
    const int64_t begin, const int64_t end)
{
  ColumnRange result;
  result.first = std::min(begin, end);
  result.second = std::max(begin, end);
  if(result.first > contig_info.m_length)
    throw GenomicsDBConfigException(std::string("Position ")+std::to_string(result.first)
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

const std::string& GenomicsDBConfigBase::get_workspace(const int rank) const
{
  VERIFY_OR_THROW((m_single_workspace_path || static_cast<size_t>(rank) < m_workspaces.size())
     && ("Workspace not defined for rank "+std::to_string(rank)).c_str());
  if(m_single_workspace_path)
    return m_workspaces[0];
  return m_workspaces[rank];
}

const std::string& GenomicsDBConfigBase::get_array_name(const int rank) const
{
  VERIFY_OR_THROW((m_single_array_name || static_cast<size_t>(rank) < m_array_names.size())
      && ("Could not find array for rank "+std::to_string(rank)).c_str());
  if(m_single_array_name)
    return m_array_names[0];
  return m_array_names[rank];
}

RowRange GenomicsDBConfigBase::get_row_partition(const int rank, const unsigned idx) const
{
  if(!m_row_partitions_specified)
    return RowRange(0, INT64_MAX-1);
  auto fixed_rank = m_single_query_row_ranges_vector ? 0 : rank;
  if(static_cast<size_t>(fixed_rank) >= m_row_ranges.size())
    throw GenomicsDBConfigException(std::string("No row partition/query interval available for process with rank ")
        +std::to_string(rank));
  VERIFY_OR_THROW(idx < m_row_ranges[fixed_rank].size());
  return m_row_ranges[fixed_rank][idx];
}

ColumnRange GenomicsDBConfigBase::get_column_partition(const int rank, const unsigned idx) const
{
  if(!m_column_partitions_specified)
    return ColumnRange(0, INT64_MAX-1);
  auto fixed_rank = m_single_query_column_ranges_vector ? 0 : rank;
  if(static_cast<size_t>(fixed_rank) >= m_column_ranges.size())
    throw GenomicsDBConfigException(std::string("No column partition/query interval available for process with rank ")
        +std::to_string(rank));
  VERIFY_OR_THROW(idx < m_column_ranges[fixed_rank].size());
  return m_column_ranges[fixed_rank][idx];
}

const std::vector<RowRange>& GenomicsDBConfigBase::get_query_row_ranges(const int rank) const
{
  auto fixed_rank = m_single_query_row_ranges_vector ? 0 : rank;
  if(static_cast<size_t>(fixed_rank) >= m_row_ranges.size())
    throw GenomicsDBConfigException(std::string("No row partition/query row range available for process with rank ")
        +std::to_string(rank));
  return m_row_ranges[fixed_rank];
}

const std::vector<ColumnRange>& GenomicsDBConfigBase::get_query_column_ranges(const int rank) const
{
  auto fixed_rank = m_single_query_column_ranges_vector ? 0 : rank;
  if(static_cast<size_t>(fixed_rank) >= m_column_ranges.size())
    throw GenomicsDBConfigException(std::string("No column partition/query column range available for process with rank ")
        +std::to_string(rank));
  return m_column_ranges[fixed_rank];
}

//Loader config functions
GenomicsDBImportConfig::GenomicsDBImportConfig()
  : GenomicsDBConfigBase()
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
  m_segment_size = 10u*1024u*1024u; //10MiB default
  m_num_cells_per_tile = 1024u;
  m_fail_if_updating = false;
  m_tiledb_compression_level = Z_DEFAULT_COMPRESSION;
  m_consolidate_tiledb_array_after_load = false;
  m_discard_missing_GTs = false;
  m_no_mandatory_VCF_fields = false;
}

void GenomicsDBConfigBase::set_vcf_header_filename(const std::string& vcf_header_filename)
{
  m_vcf_header_filename = vcf_header_filename;
  // Move m_vcf_header_filename contents to a temporary local file for non-local URIs
  if (TileDBUtils::is_cloud_path(m_vcf_header_filename)) {
    char tmp_filename[PATH_MAX];
    TileDBUtils::create_temp_filename(tmp_filename, PATH_MAX);
    TileDBUtils::move_across_filesystems(m_vcf_header_filename, tmp_filename);
    m_vcf_header_filename.assign(tmp_filename);
    m_is_tmp_vcf_header_filename = true;
  }
}

void GenomicsDBConfigBase::set_vcf_output_format(const std::string& output_format)
{
  m_vcf_output_format = output_format;
  if(m_vcf_output_format_to_is_bcf_flag.find(output_format) == m_vcf_output_format_to_is_bcf_flag.end())
  {
    std::cerr << "INFO: Invalid BCF/VCF output format: " << output_format
      <<", will output compressed VCF\n";
    m_vcf_output_format = "z";
  }
}

void GenomicsDBImportConfig::fix_callset_row_idx_bounds(const int rank)
{
  m_lb_callset_row_idx = std::max<int64_t>(m_lb_callset_row_idx, 0);
  m_ub_callset_row_idx = std::max<int64_t>(m_ub_callset_row_idx, 0);
  if(m_ub_callset_row_idx < m_lb_callset_row_idx)
    std::swap(m_lb_callset_row_idx, m_ub_callset_row_idx);
  if(m_row_based_partitioning)
  {
    m_lb_callset_row_idx = std::max(m_lb_callset_row_idx, get_row_partition(rank).first);
    m_ub_callset_row_idx = std::min(m_ub_callset_row_idx, get_row_partition(rank).second);
  }
  if(m_vid_mapper.is_initialized())
    m_ub_callset_row_idx = std::min(m_ub_callset_row_idx, m_vid_mapper.get_max_callset_row_idx());
}
 
void GenomicsDBConfigBase::update_from_loader(const GenomicsDBImportConfig& loader_config, const int rank)
{
  // Update workspace if it is not provided in GenomicsDBConfigBase
  if(m_workspaces.size() == 0u)
  {
    // Set as single workspace
    m_single_workspace_path = true;
    m_workspaces.push_back(loader_config.get_workspace(rank));
  }
  // Update array if it is not provided in query json
  if(m_array_names.size() == 0u)
  {
    // Set as single array
    m_single_array_name = true;
    m_array_names.push_back(loader_config.get_array_name(rank));
  }
  if(!m_vid_mapper.is_initialized())
  {
    m_vid_mapper = loader_config.get_vid_mapper();
    VERIFY_OR_THROW(m_vid_mapper.is_initialized());
  }
}

void GenomicsDBConfigBase::subset_query_column_ranges_based_on_partition(const GenomicsDBImportConfig& loader_config, const int rank)
{
  // Check if the partitioning is column-based, if so, pick the column corresponding to the rank
  // and update the m_column_ranges
  if (loader_config.is_partitioned_by_column())
  {
    ColumnRange my_rank_loader_column_range = loader_config.get_column_partition(rank);
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
