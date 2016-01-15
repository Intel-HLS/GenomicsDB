#ifndef VCF2BINARY_H
#define VCF2BINARY_H

#ifdef HTSDIR

#include "headers.h"
#include "vid_mapper.h"
#include "column_partition_batch.h"
#include "htslib/synced_bcf_reader.h"
#include "special_values.h"
#include "gt_common.h"

//template function for obtaining TileDB null value
template<class T>
inline T get_tiledb_null_value();

template<>
inline char get_tiledb_null_value<char>() { return NULL_CHAR; }

template<>
inline int get_tiledb_null_value<int>() { return NULL_INT; }

template<>
inline int64_t get_tiledb_null_value<int64_t>() { return NULL_INT64_T; }

template<>
inline size_t get_tiledb_null_value<size_t>() { return NULL_SIZE_T; }

template<>
inline float get_tiledb_null_value<float>() { return NULL_FLOAT; }

template<>
inline double get_tiledb_null_value<double>() { return NULL_DOUBLE; }

//Exceptions thrown 
class VCF2BinaryException : public std::exception {
  public:
    VCF2BinaryException(const std::string m="") : msg_("VCF2BinaryException : "+m) { ; }
    ~VCF2BinaryException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class VCFColumnPartition
{
  friend class VCF2Binary;
  public:
    /*
     * Primary constructor
     */
    VCFColumnPartition()
    {
      m_reader = 0;
      //buffer for vcf get functions - 16 KB
      m_vcf_get_buffer_size = 16*1024;
      m_vcf_get_buffer = (uint8_t*)malloc(m_vcf_get_buffer_size*sizeof(uint8_t));
    }
    //Delete copy constructor
    VCFColumnPartition(const VCFColumnPartition& other) = delete;
    //Define move constructor
    VCFColumnPartition(VCFColumnPartition&& other);
    ~VCFColumnPartition()
    {
      if(m_reader)
        bcf_sr_destroy(m_reader);
      m_begin_buffer_offset_for_local_callset.clear();
      m_last_full_line_end_buffer_offset_for_local_callset.clear();
      m_buffer_offset_for_local_callset.clear();
      if(m_vcf_get_buffer && m_vcf_get_buffer_size)
        free(m_vcf_get_buffer);
      m_vcf_get_buffer = 0;
      m_vcf_get_buffer_size = 0;
    }
  protected:
    int64_t m_column_interval_begin;
    int64_t m_column_interval_end;
    //Synced VCF reader
    bcf_srs_t* m_reader;
    //Position in contig from which to fetch next batch of cells
    int m_local_contig_idx;
    int64_t m_contig_position;  //position in contig (0-based)
    int64_t m_contig_tiledb_column_offset;
    //Buffer offsets - 1 per callset
    //Offset at which data should be copied for the current batch
    std::vector<int64_t> m_begin_buffer_offset_for_local_callset;
    //Offset at which the current line begins
    std::vector<int64_t> m_last_full_line_end_buffer_offset_for_local_callset;
    //Current value of offset
    std::vector<int64_t> m_buffer_offset_for_local_callset;
    //Buffer for obtaining data from htslib 
    uint8_t* m_vcf_get_buffer;
    uint64_t m_vcf_get_buffer_size;
    //Used by move constructor
    void copy_simple_members(const VCFColumnPartition& other);
};

class VCF2Binary
{
  public:
    VCF2Binary(const std::string& vcf_filename, const std::vector<std::vector<std::string>>& vcf_fields,
        unsigned file_idx, VidMapper& vid_mapper, const std::vector<ColumnRange>& partition_bounds,
        size_t max_size_per_callset,
        bool treat_deletions_as_intervals, bool parallel_partitions=false, bool noupdates=true, bool close_file=false);
    //Delete default copy constructor as it is incorrect
    VCF2Binary(const VCF2Binary& other) = delete;
    //Define move constructor explicitly
    VCF2Binary(VCF2Binary&& other);
    ~VCF2Binary();
    void clear();
    //Initialization functions
    bcf_srs_t* initialize_reader(bool open_file);
    void initialize(const std::vector<ColumnRange>& partition_bounds);
    void initialize_partition(unsigned idx, const std::vector<ColumnRange>& partition_bounds );
    /*
     * Set order of enabled callsets
     */
    void set_order_of_enabled_callsets(int64_t& order_value, std::vector<int64_t>& tiledb_row_idx_to_order) const;
    /*
     * List active row idxs
     */
    void list_active_row_idxs(const ColumnPartitionBatch& partition_batch, int64_t& row_idx_offset, std::vector<int64_t>& row_idx_vec) const;
    /*
     * Buffer
     * One entry per partition
     */
    void read_next_batch(std::vector<std::vector<uint8_t>*>& buffer_vec,
        std::vector<ColumnPartitionBatch>& partition_batches, bool close_file);
    void read_next_batch(std::vector<uint8_t>& buffer,VCFColumnPartition& vcf_partition,
        ColumnPartitionFileBatch& partition_file_batch, bool close_file);
    bool seek_and_fetch_position(VCFColumnPartition& vcf_partition, bool force_seek, bool advance_reader);
    void update_local_contig_idx(VCFColumnPartition& vcf_partition, const bcf1_t* line);
    //VCF->TileDB conversion functions
    bool convert_VCF_to_binary_for_callset(std::vector<uint8_t>& buffer, VCFColumnPartition& vcf_partition,
        size_t size_per_callset, uint64_t enabled_callsets_idx);
    /*
     * field_type_idx: 0=INFO, 1=FORMAT
     */
    template<class FieldType>
    bool convert_field_to_tiledb(std::vector<uint8_t>& buffer, VCFColumnPartition& vcf_partition, 
        int64_t& buffer_offset, const int64_t buffer_offset_limit, int local_callset_idx,
        const std::string& field_name, unsigned field_type_idx);
    /*
     * Printer
     */
    template<class FieldType>
    bool tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, const FieldType val, bool print_sep=true);
    /*
     * Null value printer
     */
    template<class FieldType>
    bool tiledb_buffer_print_null(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit);
  private:
    bool m_parallel_partitions;
    bool m_noupdates;
    bool m_close_file;
    bool m_treat_deletions_as_intervals;
    std::string m_vcf_filename;
    VidMapper* m_vid_mapper;
    int64_t m_file_idx;
    size_t m_max_size_per_callset;
    //Vector of vector of strings, outer vector has 2 elements - 0 for INFO, 1 for FORMAT
    const std::vector<std::vector<std::string>>* m_vcf_fields; 
    //Common reader when m_parallel_partitions==false
    bcf_srs_t* m_reader;
    std::string m_regions;
    //Local callset idx to tiledb row idx
    std::vector<int64_t> m_local_callset_idx_to_tiledb_row_idx;
    //Enabled local callset idx
    std::vector<int64_t> m_enabled_local_callset_idx_vec;
    //Local contig idx to global contig idx
    std::vector<int> m_local_contig_idx_to_global_contig_idx;
    //Local field idx to global field idx
    std::vector<int> m_local_field_idx_to_global_field_idx;
    std::vector<VCFColumnPartition> m_partitions;
    //Used by move constructor
    void copy_simple_members(const VCF2Binary& other);
};

#endif //ifdef HTSDIR

#endif
