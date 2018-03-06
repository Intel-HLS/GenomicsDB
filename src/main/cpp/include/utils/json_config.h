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

#ifndef RUN_CONFIG_H
#define RUN_CONFIG_H

#include "variant_query_config.h"
#include "vcf_adapter.h"
#include "vid_mapper.h"

#include "rapidjson/document.h"
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/writer.h"
#include "rapidjson/filewritestream.h"
#include "rapidjson/prettywriter.h"

//Exceptions thrown 
class RunConfigException : public std::exception {
  public:
    RunConfigException(const std::string m="") : msg_("RunConfigException : "+m) { ; }
    ~RunConfigException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

extern const char* g_json_indent_unit;

class JSONConfigBase
{
  public:
    JSONConfigBase()
    {
      m_single_array_name = false;
      m_single_workspace_path = false;
      m_single_query_column_ranges_vector = false;
      m_column_partitions_specified = false;
      m_single_query_row_ranges_vector = false;
      m_row_partitions_specified = false;
      m_scan_whole_array = false;
      //Lower and upper bounds of callset row idx to import in this invocation
      m_lb_callset_row_idx = 0;
      m_ub_callset_row_idx = INT64_MAX-1;
      clear();
    }
    void clear();
    static void extract_contig_interval_from_object(const rapidjson::Value& curr_json_object,
        const VidMapper* id_mapper, ColumnRange& result);
    static bool extract_interval_from_PB_struct_or_return_false(const rapidjson::Value& curr_json_object,
        ColumnRange& result);
    void read_from_file(const std::string& filename, const VidMapper* id_mapper=0, const int rank=0);
    const std::string& get_workspace(const int rank) const;
    const std::string& get_array_name(const int rank) const;
    ColumnRange get_column_partition(const int rank, const unsigned idx=0u) const;
    RowRange get_row_partition(const int rank, const unsigned idx=0u) const;
    const std::vector<ColumnRange> get_sorted_column_partitions() const { return m_sorted_column_partitions; }
    void read_and_initialize_vid_and_callset_mapping_if_available(FileBasedVidMapper* id_mapper, const int rank);
    const std::vector<ColumnRange>& get_query_column_ranges(const int rank) const;
    const std::vector<RowRange>& get_query_row_ranges(const int rank) const;
  protected:
    bool m_single_workspace_path;
    bool m_single_array_name;
    bool m_single_query_column_ranges_vector;
    bool m_column_partitions_specified;
    bool m_single_query_row_ranges_vector;
    bool m_row_partitions_specified;
    bool m_scan_whole_array;
    rapidjson::Document m_json;
    std::vector<std::string> m_workspaces;
    std::vector<std::string> m_array_names;
    std::vector<std::vector<ColumnRange>> m_column_ranges;
    std::vector<std::vector<RowRange>> m_row_ranges;
    std::vector<std::string> m_attributes;
    std::vector<ColumnRange> m_sorted_column_partitions;
    std::vector<RowRange> m_sorted_row_partitions;
    //Lower and upper bounds of callset row idx to import in this invocation
    int64_t m_lb_callset_row_idx;
    int64_t m_ub_callset_row_idx;
    //Vid mapping file
    std::string m_vid_mapping_file;
    //callset mapping file - if defined in upper level config file
    std::string m_callset_mapping_file;
};

class JSONLoaderConfig;

class JSONBasicQueryConfig : public JSONConfigBase
{
  public:
    JSONBasicQueryConfig() : JSONConfigBase()  { }
    void read_from_file(const std::string& filename, VariantQueryConfig& query_config, FileBasedVidMapper* id_mapper=0, int rank=0, JSONLoaderConfig* loader_config=0);
    void update_from_loader(JSONLoaderConfig* loader_config, const int rank);
    void subset_query_column_ranges_based_on_partition(const JSONLoaderConfig* loader_config, const int rank);
};

#define JSON_LOADER_PARTITION_INFO_BEGIN_FIELD_NAME "begin"
#define JSON_LOADER_PARTITION_INFO_END_FIELD_NAME "end"

class JSONLoaderConfig : public JSONConfigBase
{
  public:
    JSONLoaderConfig(bool vid_mapper_file_required = true);
    void read_from_file(const std::string& filename, FileBasedVidMapper* id_mapper=0, int rank=0);
    inline bool is_partitioned_by_row() const { return m_row_based_partitioning; }
    inline bool is_partitioned_by_column() const { return !m_row_based_partitioning; }
    inline ColumnRange get_column_partition(int idx) const
    {
      return m_row_based_partitioning ? ColumnRange(0, INT64_MAX) : JSONConfigBase::get_column_partition(idx);
    }
    inline int64_t get_max_num_rows_in_array() const { return m_max_num_rows_in_array; }
    inline bool offload_vcf_output_processing() const { return m_offload_vcf_output_processing; }
    inline bool ignore_cells_not_in_partition() const { return m_ignore_cells_not_in_partition; }
    inline bool compress_tiledb_array() const { return m_compress_tiledb_array; }
    inline bool disable_synced_writes() const { return m_disable_synced_writes; }
    inline bool delete_and_create_tiledb_array() const { return m_delete_and_create_tiledb_array; }
    inline size_t get_segment_size() const { return m_segment_size; }
    inline size_t get_num_cells_per_tile() const { return m_num_cells_per_tile; }
    inline int64_t get_tiledb_compression_level() const { return m_tiledb_compression_level; }
    inline const std::string& get_vid_mapping_filename() const { return m_vid_mapping_file; }
    inline const std::string& get_callset_mapping_filename() const { return m_callset_mapping_file; }
    inline RowRange get_row_bounds() const { return RowRange(m_lb_callset_row_idx, m_ub_callset_row_idx); }
    inline void set_vid_mapper_file_required(bool val) {
      m_vid_mapper_file_required = val;
    }
    inline bool fail_if_updating() const { return m_fail_if_updating; }
    inline bool consolidate_tiledb_array_after_load() const { return m_consolidate_tiledb_array_after_load; }
    inline bool discard_missing_GTs() const { return m_discard_missing_GTs; }
    inline bool no_mandatory_VCF_fields() const { return m_no_mandatory_VCF_fields; }
  protected:
    bool m_standalone_converter_process;
    bool m_treat_deletions_as_intervals;
    bool m_produce_combined_vcf;
    bool m_produce_tiledb_array;
    bool m_compress_tiledb_array;
    bool m_disable_synced_writes;
    bool m_delete_and_create_tiledb_array;
    bool m_row_based_partitioning;
    //do ping-pong buffering
    bool m_do_ping_pong_buffering;
    //Offload VCF output processing to another thread
    bool m_offload_vcf_output_processing;
    //Ignore cells that do not belong to this partition
    bool m_ignore_cells_not_in_partition;
    //Flag that controls whether the VCF indexes should be discarded to reduce memory consumption
    bool m_discard_vcf_index;
    unsigned m_num_entries_in_circular_buffer;
    //#VCF files to open/process in parallel
    int m_num_parallel_vcf_files;
    int m_num_converter_processes;
    int64_t m_per_partition_size;
    int64_t m_max_size_per_callset;
    //max #rows - defining domain of the array
    int64_t m_max_num_rows_in_array;
    //segment size for TileDB array
    size_t m_segment_size;
    //TileDB array #cells/tile
    size_t m_num_cells_per_tile;
    //TileDB compression level
    int m_tiledb_compression_level;
    //flag to say whether vid_mapping_file is required or optional
    bool m_vid_mapper_file_required;
    //flag that causes the loader to fail if this is an update (rather than a fresh load)
    bool m_fail_if_updating;
    //consolidate TileDB array after load - merges fragments
    bool m_consolidate_tiledb_array_after_load;
    //Discard entries with ./. or .|. as the GT field
    bool m_discard_missing_GTs;
    //The array will NOT contain mandatory VCF fields (ref, alt, qual, filter) 
    //if this flag is enabled
    bool m_no_mandatory_VCF_fields;
};

#ifdef HTSDIR

class JSONVCFAdapterConfig : public JSONConfigBase
{
  public:
    JSONVCFAdapterConfig() : JSONConfigBase()
    {
      m_vcf_header_filename = "";
      m_determine_sites_with_max_alleles = 0;
      m_combined_vcf_records_buffer_size_limit = DEFAULT_COMBINED_VCF_RECORDS_BUFFER_SIZE;
    }
    void read_from_file(const std::string& filename,
        VCFAdapter& vcf_adapter, std::string output_format="", int rank=0,
        const size_t combined_vcf_records_buffer_size_limit=0u);
    inline unsigned get_determine_sites_with_max_alleles() const { return m_determine_sites_with_max_alleles; }
    inline unsigned get_max_diploid_alt_alleles_that_can_be_genotyped() const { return m_max_diploid_alt_alleles_that_can_be_genotyped; }
    inline size_t get_combined_vcf_records_buffer_size_limit() const { return m_combined_vcf_records_buffer_size_limit; }
  protected:
    std::string m_vcf_header_filename;
    std::string m_reference_genome;
    std::string m_vcf_output_filename;
    //Count max #alt alleles , don't create combined gVCF
    unsigned m_determine_sites_with_max_alleles;
    //Max diploid alleles for which fields whose length is equal to the number of genotypes can be produced (such as PL)
    unsigned m_max_diploid_alt_alleles_that_can_be_genotyped;
    //Buffer size for combined vcf records
    size_t m_combined_vcf_records_buffer_size_limit;
};

class JSONVCFAdapterQueryConfig : public JSONVCFAdapterConfig, public JSONBasicQueryConfig
{
  public:
    JSONVCFAdapterQueryConfig() : JSONVCFAdapterConfig(), JSONBasicQueryConfig() { ; }
    void read_from_file(const std::string& filename, VariantQueryConfig& query_config,
        VCFAdapter& vcf_adapter, FileBasedVidMapper* id_mapper,
        std::string output_format="", int rank=0,
        const size_t combined_vcf_records_buffer_size_limit=0u);
};



#endif

#endif
