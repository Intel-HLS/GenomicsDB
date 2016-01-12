#ifndef TILEDB_PARTITION_LOADER_H
#define TILEDB_PARTITION_LOADER_H

#include "vid_mapper.h"
#include "column_partition_batch.h"
#include "vcf2binary.h"

//Exceptions thrown 
class VCF2TileDBException : public std::exception{
  public:
    VCF2TileDBException(const std::string m="") : msg_("VCF2TileDBException : "+m) { ; }
    ~VCF2TileDBException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

//Used to exchange info between loader and converter
//The actual buffer sizes depend on whether the object is being initialized by the loader or converter
//Using flat buffers helps MPI message passing
class LoaderConverterMessageExchange
{
  public:
    LoaderConverterMessageExchange() { m_is_serviced = false; }
    void resize_vectors(int num_divisions, int64_t total_size);
    void initialize_from_converter(int num_partitions, int64_t num_owned_callsets);
    //Used when no standalone converter processes exist
    void initialize_from_loader(int64_t all_callsets);
    inline int64_t get_idx_offset_for_converter(int converter_idx) const
    {
      assert(static_cast<size_t>(converter_idx) < m_idx_offset_per_division.size());
      return m_idx_offset_per_division[converter_idx];
    }
    inline int64_t get_idx_offset_for_partition(int partition_idx) const
    {
      assert(static_cast<size_t>(partition_idx) < m_idx_offset_per_division.size());
      return m_idx_offset_per_division[partition_idx];
    }
    bool is_partition_requested_by_loader(unsigned partition_idx) const
    { 
      assert(partition_idx < m_all_num_tiledb_row_idx_vec_request.size());
      return (m_all_num_tiledb_row_idx_vec_request[partition_idx] > 0);
    }
    bool is_new_data_in_converter_response(unsigned converter_idx) const 
    {
      assert(converter_idx < m_all_num_tiledb_row_idx_vec_response.size());
      return (m_all_num_tiledb_row_idx_vec_response[converter_idx] > 0);
    }
  public:
    bool m_is_serviced;
    //Vector containing number of tiledb row idx requested
    std::vector<int64_t> m_all_num_tiledb_row_idx_vec_request;
    //Vector containing requested row idxs
    std::vector<int64_t> m_all_tiledb_row_idx_vec_request;
    //Response by loader
    std::vector<int64_t> m_all_num_tiledb_row_idx_vec_response;
    std::vector<int64_t> m_all_tiledb_row_idx_vec_response;
  private:
    //max #values per division
    std::vector<int64_t> m_max_num_values_per_division;
    //Offset per division
    std::vector<int64_t> m_idx_offset_per_division;
};

class VCF2TileDBLoaderConverterBase
{
  public:
    VCF2TileDBLoaderConverterBase(const std::string& config_filename, int idx);
    void clear();
  protected:
    int m_idx;
    bool m_standalone_converter_process;
    bool m_treat_deletions_as_intervals;
    unsigned m_num_entries_in_circular_buffer;
    int m_num_converter_processes;
    int64_t m_per_partition_size;
    int64_t m_max_size_per_callset;
    //Vid mapping file
    std::string m_vid_mapping_filename;
    //Partition begin values
    std::vector<int64_t> m_column_partition_begin_values;
    //Ping-pong buffers
    //Note that these buffers may remain at size 0, if the ping pong buffers are owned by a different object
    std::vector<std::vector<uint8_t>> m_ping_pong_buffers;
    //Data structure for exchanging info between loader and converter
    //Note that these buffers may remain at size 0, if the exchanges are owned by a different object
    std::vector<LoaderConverterMessageExchange> m_owned_exchanges;
    void resize_circular_buffers(unsigned num_entries)
    {
      m_num_entries_in_circular_buffer = num_entries;
      m_ping_pong_buffers.resize(num_entries);
    }
};

class VCF2TileDBConverter : public VCF2TileDBLoaderConverterBase
{
  public:
    //If vid_mapper==0, build from scratch
    VCF2TileDBConverter(const std::string& config_filename, int idx,
        VidMapper* vid_mapper=0, std::vector<std::vector<uint8_t>>* buffers=0,
        std::vector<LoaderConverterMessageExchange>* exchange_vector=0);
    //Delete copy constructor
    VCF2TileDBConverter(const VCF2TileDBConverter& other) = delete;
    //Delete move constructor
    VCF2TileDBConverter(VCF2TileDBConverter&& other) = delete;
    ~VCF2TileDBConverter();
    void activate_next_batch(const unsigned exchange_idx, const int partition_idx);
    void read_next_batch(const unsigned exchange_idx);
    void dump_latest_buffer(unsigned exchange_idx, std::ostream& osptr) const;
    inline int64_t get_order_for_row_idx(const int64_t row_idx) const
    {
      assert(static_cast<size_t>(row_idx) < m_tiledb_row_idx_to_order.size());
      return m_tiledb_row_idx_to_order[row_idx];
    }
  private:
    void clear();
    void initialize_column_batch_objects();
    void initialize_vcf2binary_objects();
  private:
    int64_t m_num_callsets_owned;
    VidMapper* m_vid_mapper;
    //One per partition
    std::vector<ColumnPartitionBatch> m_partition_batch;
    //Vector of vector of strings, outer vector corresponds to FILTER, INFO, FORMAT
    std::vector<std::vector<std::string>> m_vcf_fields; 
    //One per VCF file
    std::vector<VCF2Binary> m_vcf2binary_handlers;
    //Ordering of TileDB row idx for this converter - determined by the order of m_vcf2binary_handlers
    std::vector<int64_t> m_tiledb_row_idx_to_order;
    //References to ping-pong buffers 
    //May point to buffers owned by this object or by VCF2TileDBLoader depending on the configuration
    std::vector<std::vector<uint8_t>*> m_cell_data_buffers; 
    //Data structure for exchanging info between loader and converter
    //If standalone, points to owned exchanges, else must point to those owned by VCF2TileDBLoader
    std::vector<LoaderConverterMessageExchange*> m_exchanges;
};

class CellPQElement
{
  public:
    CellPQElement()
    {
      m_offset = 0;
      m_crossed_one_buffer = false;
      m_completed = false;
    }
    bool m_crossed_one_buffer;
    bool m_completed;
    int64_t m_row_idx;
    int64_t m_column;
    int64_t m_offset;
}; 

struct TileDBCellsColumnMajorCompare
{
  bool operator()(const CellPQElement* a, const CellPQElement* b)
  {
    return ((a->m_column > b->m_column) || (a->m_column == b->m_column && a->m_row_idx > b->m_row_idx));
  }
};

typedef std::priority_queue<CellPQElement*, std::vector<CellPQElement*>, TileDBCellsColumnMajorCompare> TileDBColumnMajorPQ; 

//One per array column partition
class VCF2TileDBLoader : public VCF2TileDBLoaderConverterBase
{
  public:
    VCF2TileDBLoader(const std::string& config_filename, int idx);
    //Delete copy constructor
    VCF2TileDBLoader(const VCF2TileDBLoader& other) = delete;
    //Delete move constructor
    VCF2TileDBLoader(VCF2TileDBLoader&& other) = delete;
    ~VCF2TileDBLoader()
    {
      clear();
      if(m_converter)
        delete m_converter;
      m_converter = 0;
      delete m_vid_mapper;
      m_vid_mapper = 0;
    }
    void clear();
    VCF2TileDBConverter* get_converter() { return m_converter; }
    void read_all();
    /*
     * Get the order value at which the given row idx appears
     */
    inline int64_t get_order_for_row_idx(const int64_t row_idx) const
    { return m_standalone_converter_process ? row_idx : m_converter->get_order_for_row_idx(row_idx); }
    inline int64_t get_buffer_start_offset_for_row_idx(const int64_t row_idx) const
    { return get_order_for_row_idx(row_idx)*m_max_size_per_callset; }
    /*
     * Debug dumper
     * Return true if no more data available
     */
    bool dump_latest_buffer(unsigned exchange_idx, std::ostream& osptr);
    bool read_cell_from_buffer(const int64_t row_idx);
    bool read_next_cell_from_buffer(const int64_t row_idx);
    bool produce_cells_in_column_major_order(unsigned exchange_idx);
  private:
    void reserve_entries_in_circular_buffer(unsigned exchange_idx);
    void advance_write_idxs(unsigned exchange_idx);
    //Private members
    VidMapper* m_vid_mapper;
    //May be null
    VCF2TileDBConverter* m_converter;
    //Circular buffer logic
    std::vector<CircularBufferController> m_row_idx_to_buffer_control;
    //Vector to be used in PQ for producing cells in column major order
    std::vector<CellPQElement> m_pq_vector;
    TileDBColumnMajorPQ m_column_major_pq;
    //Row idxs not in PQ - need to be inserted in next call
    std::vector<int64_t> m_rows_not_in_pq;
};

#endif
