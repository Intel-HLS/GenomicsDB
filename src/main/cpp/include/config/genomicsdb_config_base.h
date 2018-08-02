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

#ifndef GENOMICSDB_CONFIG_BASE_H
#define GENOMICSDB_CONFIG_BASE_H

#include "vid_mapper.h"

//Exceptions thrown 
class GenomicsDBConfigException : public std::exception {
  public:
    GenomicsDBConfigException(const std::string m="") : msg_("GenomicsDBConfigException : "+m) { ; }
    ~GenomicsDBConfigException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class GenomicsDBImportConfig;

class GenomicsDBConfigBase
{
  public:
    GenomicsDBConfigBase()
    {
      m_single_array_name = false;
      m_single_workspace_path = false;
      m_single_query_column_ranges_vector = false;
      m_column_partitions_specified = false;
      m_single_query_row_ranges_vector = false;
      m_row_partitions_specified = false;
      m_scan_whole_array = false;
      m_produce_GT_field = false;
      m_produce_FILTER_field = false;
      m_index_output_VCF = false;
      m_sites_only_query = false;
      m_produce_GT_with_min_PL_value_for_spanning_deletions = false;
      //Lower and upper bounds of callset row idx to import in this invocation
      m_lb_callset_row_idx = 0;
      m_ub_callset_row_idx = INT64_MAX-1;
      m_segment_size = 10u*1024u*1024u; //10MiB default
      m_disable_file_locking_in_tiledb = false;
      m_determine_sites_with_max_alleles = false;
      m_max_diploid_alt_alleles_that_can_be_genotyped = MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED;
      m_combined_vcf_records_buffer_size_limit = 10*1024u;
    }
    const std::string& get_workspace(const int rank) const;
    const std::string& get_array_name(const int rank) const;
    ColumnRange get_column_partition(const int rank, const unsigned idx=0u) const;
    RowRange get_row_partition(const int rank, const unsigned idx=0u) const;
    const std::vector<ColumnRange> get_sorted_column_partitions() const { return m_sorted_column_partitions; }
    const std::vector<ColumnRange>& get_query_column_ranges(const int rank) const;
    const std::vector<RowRange>& get_query_row_ranges(const int rank) const;
    inline size_t get_segment_size() const { return m_segment_size; }
    void set_segment_size(const size_t v) { m_segment_size = v; }
    inline unsigned get_determine_sites_with_max_alleles() const { return m_determine_sites_with_max_alleles; }
    inline unsigned get_max_diploid_alt_alleles_that_can_be_genotyped() const { return m_max_diploid_alt_alleles_that_can_be_genotyped; }
    void set_combined_vcf_records_buffer_size_limit(const size_t val) { m_combined_vcf_records_buffer_size_limit = val; }
    inline size_t get_combined_vcf_records_buffer_size_limit() const { return m_combined_vcf_records_buffer_size_limit; }
    void set_vcf_header_filename(const std::string& vcf_header_filename);
    const std::string& get_vcf_header_filename() const { return m_vcf_header_filename; }
    void set_vcf_output_format(const std::string& output_format);
    const std::string& get_vcf_output_format() const { return m_vcf_output_format; }
    const std::string& get_vcf_output_filename() const { return m_vcf_output_filename; }
    const std::string& get_reference_genome() const { return m_reference_genome; }
    const bool produce_GT_field() const { return m_produce_GT_field; }
    const bool produce_FILTER_field() const { return m_produce_FILTER_field; }
    const bool sites_only_query() const { return m_sites_only_query; }
    const bool index_output_VCF() const { return m_index_output_VCF; }
    const bool produce_GT_with_min_PL_value_for_spanning_deletions() const
    { return m_produce_GT_with_min_PL_value_for_spanning_deletions; }
    const VidMapper& get_vid_mapper() const { return m_vid_mapper; }
    //Utility functions
    static ColumnRange verify_contig_position_and_get_tiledb_column_interval(const ContigInfo& contig_info,
        const int64_t begin, const int64_t end);
    const std::string& get_callset_mapping_file() const { return m_callset_mapping_file; }
    const std::string& get_vid_mapping_file() const { return m_vid_mapping_file; }
    //Sometimes information is present in the loader - copy over
    void update_from_loader(const GenomicsDBImportConfig& loader_config, const int rank);
    void subset_query_column_ranges_based_on_partition(const GenomicsDBImportConfig& loader_config, const int rank);
    inline RowRange get_row_bounds() const { return RowRange(m_lb_callset_row_idx, m_ub_callset_row_idx); }
    inline uint64_t get_num_rows_within_bounds() const { return m_ub_callset_row_idx - m_lb_callset_row_idx + 1ull; }
    inline bool disable_file_locking_in_tiledb() const { return m_disable_file_locking_in_tiledb; }
  protected:
    bool m_single_workspace_path;
    bool m_single_array_name;
    bool m_single_query_column_ranges_vector;
    bool m_column_partitions_specified;
    bool m_single_query_row_ranges_vector;
    bool m_row_partitions_specified;
    bool m_scan_whole_array;
    //GATK CombineGVCF does not produce GT field by default - option to produce GT
    bool m_produce_GT_field;
    //GATK CombineGVCF does not produce FILTER field by default - option to produce FILTER
    bool m_produce_FILTER_field;
    //index output VCF file
    bool m_index_output_VCF;
    //sites-only query - doesn't produce any of the FORMAT fields
    bool m_sites_only_query;
    //when producing GT, use the min PL value GT for spanning deletions
    bool m_produce_GT_with_min_PL_value_for_spanning_deletions;
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
    //TileDB segment size
    size_t m_segment_size;
    //VCF output parameters 
    std::string m_vcf_header_filename;
    std::string m_reference_genome;
    std::string m_vcf_output_filename;
    std::string m_vcf_output_format;
    //Count max #alt alleles , don't create combined gVCF
    unsigned m_determine_sites_with_max_alleles;
    //Max diploid alleles for which fields whose length is equal to the number of genotypes can be produced (such as PL)
    unsigned m_max_diploid_alt_alleles_that_can_be_genotyped;
    //Buffer size for combined vcf records
    size_t m_combined_vcf_records_buffer_size_limit;
    //VidMapper
    VidMapper m_vid_mapper;
    //Might be empty strings if using Protobuf
    std::string m_vid_mapping_file;
    std::string m_callset_mapping_file;
    //Disable file locking in TileDB
    bool m_disable_file_locking_in_tiledb;
  public:
    //Static convenience member
    static std::unordered_map<std::string, bool> m_vcf_output_format_to_is_bcf_flag;
};

class GenomicsDBImportConfig : public GenomicsDBConfigBase
{
  public:
    GenomicsDBImportConfig();
    void read_from_file(const std::string& filename, int rank=0);
    inline bool is_partitioned_by_row() const { return m_row_based_partitioning; }
    inline bool is_partitioned_by_column() const { return !m_row_based_partitioning; }
    inline int64_t get_max_num_rows_in_array() const { return m_max_num_rows_in_array; }
    inline bool offload_vcf_output_processing() const { return m_offload_vcf_output_processing; }
    inline bool ignore_cells_not_in_partition() const { return m_ignore_cells_not_in_partition; }
    inline bool compress_tiledb_array() const { return m_compress_tiledb_array; }
    inline bool disable_synced_writes() const { return m_disable_synced_writes; }
    inline bool delete_and_create_tiledb_array() const { return m_delete_and_create_tiledb_array; }
    inline size_t get_segment_size() const { return m_segment_size; }
    inline size_t get_num_cells_per_tile() const { return m_num_cells_per_tile; }
    inline int64_t get_tiledb_compression_level() const { return m_tiledb_compression_level; }
    inline bool fail_if_updating() const { return m_fail_if_updating; }
    inline bool consolidate_tiledb_array_after_load() const { return m_consolidate_tiledb_array_after_load; }
    inline bool discard_missing_GTs() const { return m_discard_missing_GTs; }
    inline bool no_mandatory_VCF_fields() const { return m_no_mandatory_VCF_fields; }
    inline bool treat_deletions_as_intervals() const { return m_treat_deletions_as_intervals; }
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
    //flag that causes the loader to fail if this is an update (rather than a fresh load)
    bool m_fail_if_updating;
    //consolidate TileDB array after load - merges fragments
    bool m_consolidate_tiledb_array_after_load;
    //Discard entries with ./. or .|. as the GT field
    bool m_discard_missing_GTs;
    //The array will NOT contain mandatory VCF fields (ref, alt, qual, filter) 
    //if this flag is enabled
    bool m_no_mandatory_VCF_fields;
  protected:
    void fix_callset_row_idx_bounds(const int rank);
};

#endif
