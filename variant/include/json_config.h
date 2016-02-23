#ifndef RUN_CONFIG_H
#define RUN_CONFIG_H

#include "variant_query_config.h"
#include "vcf_adapter.h"
#include "vid_mapper.h"

#include "rapidjson/document.h"
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"

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
      clear();
    }
    void clear();
    void read_from_file(const std::string& filename);
    const std::string& get_workspace(const int rank) const;
    const std::string& get_array_name(const int rank) const;
    ColumnRange get_column_partition(const int rank, const unsigned idx=0u) const;
    RowRange get_row_partition(const int rank, const unsigned idx=0u) const;
    const std::vector<ColumnRange> get_sorted_column_partitions() const { return m_sorted_column_partitions; }
    std::pair<std::string, std::string> get_vid_mapping_filename(FileBasedVidMapper* id_mapper, const int rank);
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
};

class JSONBasicQueryConfig : public JSONConfigBase
{
  public:
    JSONBasicQueryConfig() : JSONConfigBase()  { }
    void read_from_file(const std::string& filename, VariantQueryConfig& query_config, FileBasedVidMapper* id_mapper=0, int rank=0);
};

class JSONLoaderConfig : public JSONConfigBase
{
  public:
    JSONLoaderConfig();
    void read_from_file(const std::string& filename, FileBasedVidMapper* id_mapper=0, int rank=0);
    inline bool is_partitioned_by_row() const { return m_row_based_partitioning; }
    inline bool is_partitioned_by_column() const { return !m_row_based_partitioning; }
    inline ColumnRange get_column_partition(int idx) const
    {
      return m_row_based_partitioning ? ColumnRange(0, INT64_MAX) : JSONConfigBase::get_column_partition(idx);
    }
  protected:
    bool m_standalone_converter_process;
    bool m_treat_deletions_as_intervals;
    bool m_produce_combined_vcf;
    bool m_produce_tiledb_array;
    bool m_row_based_partitioning;
    //Flag that controls whether the VCF indexes should be discarded to reduce memory consumption
    bool m_discard_vcf_index;
    unsigned m_num_entries_in_circular_buffer;
    int m_num_converter_processes;
    int64_t m_per_partition_size;
    int64_t m_max_size_per_callset;
    //Vid mapping file
    std::string m_vid_mapping_filename;
    //callset mapping file - if defined in upper level config file
    std::string m_callset_mapping_file;
    //Limit callset row idx to this value
    int64_t m_limit_callset_row_idx;
    //#VCF files to open/process in parallel
    int m_num_parallel_vcf_files;
    //do ping-pong buffering
    bool m_do_ping_pong_buffering;
    //Offload VCF output processing to another thread
    bool m_offload_vcf_output_processing;
};

#ifdef HTSDIR

class JSONVCFAdapterConfig : public JSONConfigBase
{
  public:
    JSONVCFAdapterConfig() : JSONConfigBase()
    {
      m_vcf_header_filename = "";
    }
    void read_from_file(const std::string& filename,
        VCFAdapter& vcf_adapter, std::string output_format="", int rank=0);
  protected:
    std::string m_vcf_header_filename;
    std::string m_reference_genome;
    std::string m_vcf_output_filename;
};

class JSONVCFAdapterQueryConfig : public JSONVCFAdapterConfig, public JSONBasicQueryConfig
{
  public:
    JSONVCFAdapterQueryConfig() : JSONVCFAdapterConfig(), JSONBasicQueryConfig() { ; }
    void read_from_file(const std::string& filename, VariantQueryConfig& query_config,
        VCFAdapter& vcf_adapter, FileBasedVidMapper* id_mapper,
        std::string output_format="", int rank=0);
};



#endif

#endif
