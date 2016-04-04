#ifndef TILEDB_IMPORTER_H
#define TILEDB_IMPORTER_H

#include "gt_common.h"
#include "column_partition_batch.h"
#include "vid_mapper.h"
#include "histogram.h"

#define PRODUCE_BINARY_CELLS 1
/*#define PRODUCE_CSV_CELLS 1*/

//Binary gets priority
#if defined(PRODUCE_BINARY_CELLS) and defined(PRODUCE_CSV_CELLS)
#undef PRODUCE_CSV_CELLS
#endif

//If none defined, produce binary cells
#if !defined(PRODUCE_BINARY_CELLS) and !defined(PRODUCE_CSV_CELLS)
#define PRODUCE_BINARY_CELLS 1
#endif

//Exceptions thrown 
class File2TileDBBinaryException : public std::exception {
  public:
    File2TileDBBinaryException(const std::string m="") : msg_("File2TileDBBinaryException : "+m) { ; }
    ~File2TileDBBinaryException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class FileReaderBase
{
  public:
    FileReaderBase() { m_is_record_valid = false; }
    virtual ~FileReaderBase() { }
    virtual void add_reader() = 0;
    virtual void remove_reader() = 0;
    virtual void read_and_advance() = 0;
  protected:
    bool m_is_record_valid;
};

class File2TileDBBinaryColumnPartitionBase
{
  friend class File2TileDBBinaryBase;
  public:
    File2TileDBBinaryColumnPartitionBase()
    {
      m_base_reader_ptr = 0;
      m_buffer_ptr = 0;
    }
    virtual ~File2TileDBBinaryColumnPartitionBase();
    void clear()
    {
      m_begin_buffer_offset_for_local_callset.clear();
      m_last_full_line_end_buffer_offset_for_local_callset.clear();
      m_buffer_offset_for_local_callset.clear();
      m_buffer_full_for_local_callset.clear();
    }
    //Delete copy constructor
    File2TileDBBinaryColumnPartitionBase(const File2TileDBBinaryColumnPartitionBase& other) = delete;
    //Define move constructor
    File2TileDBBinaryColumnPartitionBase(File2TileDBBinaryColumnPartitionBase&& other);
    void initialize_base_class_members(const int64_t begin, const int64_t end,
        const uint64_t num_enabled_callsets, FileReaderBase* ptr);
    FileReaderBase* get_base_reader_ptr() { return m_base_reader_ptr; }
    /*
     * Buffer control access functions
     */
    bool is_buffer_full(unsigned idx_in_file) const
    {
      assert(idx_in_file < m_buffer_full_for_local_callset.size());
      return m_buffer_full_for_local_callset[idx_in_file];
    }
    void set_buffer_full_if_true(unsigned idx_in_file, bool val)
    {
      assert(idx_in_file < m_buffer_full_for_local_callset.size());
      m_buffer_full_for_local_callset[idx_in_file] = m_buffer_full_for_local_callset[idx_in_file] || val;
    }
    void reset_buffer_full(unsigned idx_in_file)
    {
      assert(idx_in_file < m_buffer_full_for_local_callset.size());
      m_buffer_full_for_local_callset[idx_in_file] = false;
    }
    std::vector<uint8_t>& get_buffer()
    {
      assert(m_buffer_ptr);
      return *m_buffer_ptr;
    }
    /*
     * abstract virtual functions
     */
    virtual int64_t get_column_position_in_record() const = 0;
  protected:
    int64_t m_column_interval_begin;
    int64_t m_column_interval_end;
    //Buffer offsets - 1 per callset
    //Offset at which data should be copied for the current batch
    std::vector<int64_t> m_begin_buffer_offset_for_local_callset;
    //Offset at which the current line begins
    std::vector<int64_t> m_last_full_line_end_buffer_offset_for_local_callset;
    //Current value of offset
    std::vector<int64_t> m_buffer_offset_for_local_callset;
    //Buffer full flags - 1 per callset
    std::vector<bool> m_buffer_full_for_local_callset;
    //Pointer to buffer
    std::vector<uint8_t>* m_buffer_ptr;
    FileReaderBase* m_base_reader_ptr;
};

/*
 * Base class for all instances of objects that convert file formats for importing
 * to TileDB
 */
class File2TileDBBinaryBase
{
  public:
    File2TileDBBinaryBase(const std::string& filename,
        unsigned file_idx, VidMapper& vid_mapper,
        size_t max_size_per_callset,
        bool treat_deletions_as_intervals,
        bool parallel_partitions=false, bool noupdates=true, bool close_file=false);
    void initialize_base_column_partitions(const std::vector<ColumnRange>& partition_bounds);
    //Delete copy constructor
    File2TileDBBinaryBase(const File2TileDBBinaryBase& other) = delete;
    //Define move constructor
    File2TileDBBinaryBase(File2TileDBBinaryBase&& other);
    //Destructor
    virtual ~File2TileDBBinaryBase();
    void clear();
    /*
     * Set order of enabled callsets
     */
    void set_order_of_enabled_callsets(int64_t& order_value, std::vector<int64_t>& tiledb_row_idx_to_order) const;
    /*
     * List active row idxs
     */
    void list_active_row_idxs(const ColumnPartitionBatch& partition_batch, int64_t& row_idx_offset, std::vector<int64_t>& row_idx_vec) const;
    /*
     * */
    void read_next_batch(std::vector<std::vector<uint8_t>*>& buffer_vec,
        std::vector<ColumnPartitionBatch>& partition_batches, bool close_file);
    void read_next_batch(std::vector<uint8_t>& buffer,File2TileDBBinaryColumnPartitionBase& partition_info,
        ColumnPartitionFileBatch& partition_file_batch, bool close_file);
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
    /*
     * Create histogram
     */
    void create_histogram(uint64_t max_histogram_range, unsigned num_bins);
    UniformHistogram* get_histogram() { return m_histogram; }
    //Functions that must be over-ridden by all sub-classes
    /*
     * Initialization of column partitions by sub class
     */
    virtual void initialize_column_partitions(const std::vector<ColumnRange>& partition_bounds) = 0;
    /*
     * Create the subclass of File2TileDBBinaryColumnPartitionBase that must be used
     */
    virtual File2TileDBBinaryColumnPartitionBase* create_new_column_partition_object() const = 0;
    /*
     * Create the subclass of FileReaderBase that must be used
     */
    virtual FileReaderBase* create_new_reader_object(const std::string& filename, bool open_file) const = 0;
    /*
     * Convert current record to TileDB binary in the buffer
     */
    virtual bool convert_record_to_binary(std::vector<uint8_t>& buffer,
        File2TileDBBinaryColumnPartitionBase& partition_info) = 0;
    /*
     * Seek and/or advance to position in the file as described by partition_info
     */
    virtual bool seek_and_fetch_position(File2TileDBBinaryColumnPartitionBase& partition_info, bool force_seek,
        bool advance_reader) = 0;
    /*
     * Return #callsets in current record
     */
    virtual uint64_t get_num_callsets_in_record(const File2TileDBBinaryColumnPartitionBase& partition_info) const = 0;
  protected:
    inline int64_t get_enabled_idx_for_local_callset_idx(int64_t local_callset_idx) const
    {
      assert(local_callset_idx >= 0 && static_cast<size_t>(local_callset_idx) < m_local_callset_idx_to_enabled_idx.size());
      return m_local_callset_idx_to_enabled_idx[local_callset_idx];
    }
  protected:
    bool m_parallel_partitions;
    bool m_noupdates;
    bool m_close_file;
    bool m_treat_deletions_as_intervals;
    VidMapper* m_vid_mapper;
    int64_t m_file_idx;
    size_t m_max_size_per_callset;
    std::string m_filename;
    //Local callset idx to tiledb row idx
    std::vector<int64_t> m_local_callset_idx_to_tiledb_row_idx;
    //Enabled local callset idx
    std::vector<int64_t> m_enabled_local_callset_idx_vec;
    //Local callset idx to enabled idx
    std::vector<int64_t> m_local_callset_idx_to_enabled_idx;
    //Reader
    FileReaderBase* m_base_reader_ptr;
    //Partition read state - pointers to objects of sub-classes of File2TileDBBinaryColumnPartitionBase
    //Must be initialized by sub-classes of File2TileDBBinaryBase
    std::vector<File2TileDBBinaryColumnPartitionBase*> m_base_partition_ptrs;
    //Histogram
    UniformHistogram* m_histogram;
  private:
    //Called by move constructor
    void copy_simple_members(const File2TileDBBinaryBase& other);
};

#endif
