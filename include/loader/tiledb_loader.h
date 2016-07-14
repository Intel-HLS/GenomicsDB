/**
 * The MIT License (MIT)
 * Copyright (c) 2016 Intel Corporation
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

#ifndef TILEDB_PARTITION_LOADER_H
#define TILEDB_PARTITION_LOADER_H

#include "vid_mapper.h"
#include "column_partition_batch.h"
#include "tiledb_loader_file_base.h"
#include "load_operators.h"
#include "json_config.h"

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

class VCF2TileDBLoaderConverterBase : public JSONLoaderConfig
{
  public:
    VCF2TileDBLoaderConverterBase(const std::string& config_filename, int idx);
    inline int64_t get_column_partition_end() const
    {
      return JSONLoaderConfig::get_column_partition(m_idx).second;
    }
    inline int64_t get_column_partition_begin() const
    {
      return JSONLoaderConfig::get_column_partition(m_idx).first;
    }
    inline ColumnRange get_column_partition() const { return JSONLoaderConfig::get_column_partition(m_idx); }
    void clear();
  protected:
    void resize_circular_buffers(unsigned num_entries)
    {
      m_num_entries_in_circular_buffer = num_entries;
      m_ping_pong_buffers.resize(num_entries);
    }
    void determine_num_callsets_owned(const VidMapper* vid_mapper, const bool from_loader);
  protected:
    int m_idx;
    //Ping-pong buffers
    //Note that these buffers may remain at size 0, if the ping pong buffers are owned by a different object
    std::vector<std::vector<uint8_t>> m_ping_pong_buffers;
    //Data structure for exchanging info between loader and converter
    //Note that these buffers may remain at size 0, if the exchanges are owned by a different object
    std::vector<LoaderConverterMessageExchange> m_owned_exchanges;
    //#callsets owned by this entity
    int64_t m_num_callsets_owned;
    //One entry for every owned file - contains #callsets per owned file
    std::vector<int64_t> m_num_callsets_in_owned_file;
    //Row idxs for callsets owned by this entity - may not be initialized
    std::vector<int64_t> m_owned_row_idx_vec;
};

#ifdef HTSDIR
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
    inline int64_t get_designated_row_idx_for_order(const int64_t order) const
    {
      assert(order >= 0 && static_cast<size_t>(order) < m_order_to_designated_tiledb_row_idx.size());
      auto val = m_order_to_designated_tiledb_row_idx[order];
      assert(val >= 0 && static_cast<size_t>(val) < m_tiledb_row_idx_to_order.size());
      return val;
    }
    inline size_t get_num_order_values() const { return m_order_to_designated_tiledb_row_idx.size(); }
    void create_and_print_histogram(const std::string& config_filename, std::ostream& fptr=std::cout);
  private:
    void clear();
    void initialize_column_batch_objects();
    void initialize_file2binary_objects();
    File2TileDBBinaryBase* create_file2tiledb_object(const FileInfo& file_info, const uint64_t local_file_idx,
        const std::vector<ColumnRange>& partition_bounds);
  private:
    VidMapper* m_vid_mapper;
    //One per partition
    std::vector<ColumnPartitionBatch> m_partition_batch;
    //Vector of vector of strings, outer vector corresponds to FILTER, INFO, FORMAT
    std::vector<std::vector<std::string>> m_vcf_fields;
    //One per VCF file
    std::vector<File2TileDBBinaryBase*> m_file2binary_handlers;
    //Multiple row_idx could point to the same order
    //Ordering of TileDB row idx for this converter - determined by the order of m_file2binary_handlers
    std::vector<int64_t> m_tiledb_row_idx_to_order;
    //Order to row_idx - points to one designated row_idx for a given order value
    std::vector<int64_t> m_order_to_designated_tiledb_row_idx;
    //References to ping-pong buffers 
    //May point to buffers owned by this object or by VCF2TileDBLoader depending on the configuration
    std::vector<std::vector<uint8_t>*> m_cell_data_buffers; 
    //Data structure for exchanging info between loader and converter
    //If standalone, points to owned exchanges, else must point to those owned by VCF2TileDBLoader
    std::vector<LoaderConverterMessageExchange*> m_exchanges;
};
#endif

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
      for(auto op : m_operators)
        if(op)
          delete op;
      clear();
#ifdef HTSDIR
      if(m_converter)
        delete m_converter;
      m_converter = 0;
#endif
      delete m_vid_mapper;
      m_vid_mapper = 0;
    }
    void clear();
#ifdef HTSDIR
    VCF2TileDBConverter* get_converter() { return m_converter; }
    /*
     * For non-standalone converter processes, this function will do the read, convert and processing
     * of all VCF files
     */
    void read_all();
#endif
    /*
     * Get the order value at which the given row idx appears
     */
    inline int64_t get_order_for_row_idx(const int64_t row_idx) const
    {
#ifdef HTSDIR
      return m_standalone_converter_process ? row_idx : m_converter->get_order_for_row_idx(row_idx);
#else
      return row_idx;
#endif
    }
    inline int64_t get_designated_row_idx_for_order(const int64_t order) const
    {
#ifdef HTSDIR
      return m_standalone_converter_process ? order : m_converter->get_designated_row_idx_for_order(order);
#else
      return order;
#endif
    }
    inline size_t get_num_order_values() const
    {
#ifdef HTSDIR
      return m_standalone_converter_process ? m_num_callsets_owned : m_converter->get_num_order_values();
#else
      return m_num_callsets_owned;
#endif
    }
    inline int64_t get_buffer_start_offset_for_row_idx(const int64_t row_idx) const
    {
      assert(get_order_for_row_idx(row_idx) >= 0 && get_order_for_row_idx(row_idx) < m_num_callsets_owned);
      return get_order_for_row_idx(row_idx)*m_max_size_per_callset;
    }
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
#ifdef HTSDIR
    //May be null
    VCF2TileDBConverter* m_converter;
#endif
    //Circular buffer logic
    std::vector<CircularBufferController> m_order_idx_to_buffer_control;
    //Vector to be used in PQ for producing cells in column major order
    std::vector<CellPQElement> m_pq_vector;
    TileDBColumnMajorPQ m_column_major_pq;
    //Row idxs not in PQ - need to be inserted in next call
    std::vector<int64_t> m_designated_rows_not_in_pq;
    //Operators - act on one cell per call
    std::vector<LoaderOperatorBase*> m_operators;
    std::vector<bool> m_operators_overflow;
    unsigned m_num_operators_overflow_in_last_round;
};

#endif
