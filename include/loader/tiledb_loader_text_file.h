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

/*
 * We use libcsv (https://sourceforge.net/projects/libcsv/) to parse CSV files. libcsv
 * is licensed under the GNU Library or Lesser General Public License version 2.0 (LGPLv2).
 * So, if you are re-distributing binaries or object files, they may be subject to LGPLv2 terms.
 * Please ensure that any binaries/object files you distribute are compliant with LGPLv2.
 */ 

#ifndef TILEDB_LOADER_TEXT_FILE_H
#define TILEDB_LOADER_TEXT_FILE_H

#include "tiledb_loader_file_base.h"

//Use libcsv for parsing
#if (defined LIBCSV_DIR) && !(defined USE_LIBCSV)
#define USE_LIBCSV 1
#endif

#ifdef USE_LIBCSV
#include <csv.h>
#endif

//Exceptions thrown
class LineBasedTextFileException : public std::exception {
  public:
    LineBasedTextFileException(const std::string m="") : msg_("LineBasedTextFileException : "+m) { ; }
    ~LineBasedTextFileException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

//Line based reader
class LineBasedTextFileReader : public FileReaderBase
{
  public:
    LineBasedTextFileReader();
    //Delete copy constructor
    LineBasedTextFileReader(const LineBasedTextFileReader& other) = delete;
    //Delete move constructor
    LineBasedTextFileReader(LineBasedTextFileReader&& other) = delete;
    ~LineBasedTextFileReader();
    void initialize(const char* filename, bool open_file);
    void add_reader();
    void remove_reader();
    void read_and_advance();
    void seek(const fpos_t& pos)
    {
      assert(m_fptr);
      auto status = fsetpos(m_fptr, &pos);
      assert(status == 0);
    }
    void get_position(fpos_t& pos) const
    {
      assert(m_fptr);
      auto status = fgetpos(m_fptr, &pos);
      assert(status == 0);
    }
    inline const char* get_line() const { return m_is_record_valid ? m_line_buffer : 0; }
    inline size_t get_line_length() const  { return m_is_record_valid ? m_line_length : 0ull; }
  private:
    FILE* m_fptr;
    char* m_line_buffer;
    size_t m_line_buffer_size;
    //Line length excluding newline
    size_t m_line_length;
};

/*
 * Abstract base class for text formats - only contains fpos_t
 */
class LineBasedTextFile2TileDBBinaryColumnPartition : public File2TileDBBinaryColumnPartitionBase
{
  public:
    LineBasedTextFile2TileDBBinaryColumnPartition() : File2TileDBBinaryColumnPartitionBase()
    {
      m_initialized_file_position_to_partition_begin = false;
      memset(&m_file_position, 0, sizeof(fpos_t));
    }
    //Delete copy constructor
    LineBasedTextFile2TileDBBinaryColumnPartition(const LineBasedTextFile2TileDBBinaryColumnPartition& other) = delete;
    //Define move constructor
    LineBasedTextFile2TileDBBinaryColumnPartition(LineBasedTextFile2TileDBBinaryColumnPartition&& other)
      : File2TileDBBinaryColumnPartitionBase(std::move(other))
    {
      m_initialized_file_position_to_partition_begin = other.m_initialized_file_position_to_partition_begin;
      m_file_position = other.m_file_position;
    }
    virtual ~LineBasedTextFile2TileDBBinaryColumnPartition() { }
    bool is_initialized_file_position_to_partition_begin() const { return m_initialized_file_position_to_partition_begin; }
    void set_initialized_file_position_to_partition_begin(const bool val)
    { m_initialized_file_position_to_partition_begin = val; }
  protected:
    bool m_initialized_file_position_to_partition_begin;
    fpos_t m_file_position;
};

class CSV2TileDBBinaryColumnPartition : public LineBasedTextFile2TileDBBinaryColumnPartition
{
  friend class CSV2TileDBBinary;
  public:
    CSV2TileDBBinaryColumnPartition() : LineBasedTextFile2TileDBBinaryColumnPartition()
    {
      m_current_column_position = -1;
#ifdef USE_LIBCSV
      csv_init(&m_csv_parser, CSV_STRICT|CSV_APPEND_NULL);
#else
      throw LineBasedTextFileException("Cannot import CSV files without libcsv - recompile GenomicsDB with USE_LIBCSV and/or LIBCSV_DIR set");
#endif
    }
    //Delete copy constructor
    CSV2TileDBBinaryColumnPartition(const CSV2TileDBBinaryColumnPartition& other) = delete;
    //Define move constructor
    CSV2TileDBBinaryColumnPartition(CSV2TileDBBinaryColumnPartition&& other)
      : LineBasedTextFile2TileDBBinaryColumnPartition(std::move(other))
    {
      m_current_column_position = other.m_current_column_position;
#ifdef USE_LIBCSV
      std::swap(m_csv_parser, other.m_csv_parser);
#endif
    }
    ~CSV2TileDBBinaryColumnPartition()
    {
#ifdef USE_LIBCSV
      csv_free(&m_csv_parser);
#endif
    }
    int64_t get_column_position_in_record() const { return m_current_column_position; }
  private:
    int64_t m_current_column_position;
#ifdef USE_LIBCSV
    struct csv_parser m_csv_parser;
#endif
};

/*
 * Abstract base class for all instances of text file converters
 */
class LineBasedTextFile2TileDBBinary : public File2TileDBBinaryBase
{
  public:
    LineBasedTextFile2TileDBBinary(const std::string& filename,
        unsigned file_idx, VidMapper& vid_mapper,
        size_t max_size_per_callset,
        bool treat_deletions_as_intervals,
        bool parallel_partitions=false, bool noupdates=true, bool close_file=false)
      : File2TileDBBinaryBase(filename, file_idx, vid_mapper,
          max_size_per_callset,
          treat_deletions_as_intervals,
          parallel_partitions, noupdates, close_file)
    {
      vid_mapper.build_tiledb_array_schema(m_array_schema, "dummy", false, RowRange(0, INT64_MAX-1), false);
    }
    //Delete copy constructor
    LineBasedTextFile2TileDBBinary(const LineBasedTextFile2TileDBBinary& other) = delete;
    //Define move constructor
    LineBasedTextFile2TileDBBinary(LineBasedTextFile2TileDBBinary&& other)
      : File2TileDBBinaryBase(std::move(other))
    {
      std::swap(m_array_schema, other.m_array_schema);
    }
    virtual ~LineBasedTextFile2TileDBBinary()
    {
      if(m_array_schema)
        delete m_array_schema;
      m_array_schema = 0;
    }
    /*
     * Create the subclass of GenomicsDBImportReaderBase that must be used
     */
    GenomicsDBImportReaderBase* create_new_reader_object(const std::string& filename, bool open_file) const
    {
      return dynamic_cast<GenomicsDBImportReaderBase*>(new LineBasedTextFileReader());
    }
  protected:
    VariantArraySchema* m_array_schema;
};

enum TileDBCSVFieldPosIdxEnum
{
  TILEDB_CSV_ROW_POS_IDX=0u,
  TILEDB_CSV_COLUMN_POS_IDX,
  TILEDB_CSV_END_IDX,
};

//Forward declaration for pointer
class CSVLineParseStruct;

//Callback functions
void csv_parse_callback(void* field_ptr, size_t field_size, void* parse_ptr);
void csv_line_end_callback(int terminating_token, void* parse_ptr);

class CSV2TileDBBinary : public LineBasedTextFile2TileDBBinary
{
  public:
    CSV2TileDBBinary(const std::string& filename,
        unsigned file_idx, VidMapper& vid_mapper,
        size_t max_size_per_callset,
        const std::vector<ColumnRange>& partition_bounds,
        bool treat_deletions_as_intervals,
        bool parallel_partitions=false, bool noupdates=true, bool close_file=false);
    //Delete copy constructor
    CSV2TileDBBinary(const CSV2TileDBBinary& other) = delete;
    //Define move constructor
    CSV2TileDBBinary(CSV2TileDBBinary&& other)
      : LineBasedTextFile2TileDBBinary(std::move(other))
    {
      m_cleanup_file = other.m_cleanup_file;
      other.m_cleanup_file = false;
    }
    ~CSV2TileDBBinary();
    //Functions that must be over-ridden by all sub-classes
    /*
     * Set order of enabled callsets
     * In a CSV, a line contains information from a single callset
     */
    void set_order_of_enabled_callsets(int64_t& order_value, std::vector<int64_t>& tiledb_row_idx_to_order) const;
    /*
     * List active row idxs
     * In a CSV, a line contains information from a single callset
     */
    void list_active_row_idxs(const ColumnPartitionBatch& partition_batch, int64_t& row_idx_offset, std::vector<int64_t>& row_idx_vec) const;
    /*
     * Initialization of column partitions by sub class
     */
    void initialize_column_partitions(const std::vector<ColumnRange>& partition_bounds);
    /*
     * Create the subclass of File2TileDBBinaryColumnPartitionBase that must be used
     */
    File2TileDBBinaryColumnPartitionBase* create_new_column_partition_object() const
    {
      return dynamic_cast<File2TileDBBinaryColumnPartitionBase*>(new CSV2TileDBBinaryColumnPartition());
    }
    /*
     * Convert current record to TileDB binary in the buffer
     */
    bool convert_record_to_binary(std::vector<uint8_t>& buffer,
        File2TileDBBinaryColumnPartitionBase& partition_info)
    {
      auto& csv_partition_info = dynamic_cast<CSV2TileDBBinaryColumnPartition&>(partition_info);
      auto csv_reader_ptr = dynamic_cast<LineBasedTextFileReader*>(partition_info.get_base_reader_ptr());
      assert(csv_reader_ptr);
      assert(csv_reader_ptr->get_line());
      return parse_line(csv_reader_ptr->get_line(), csv_partition_info, UINT32_MAX, true);
    }
    /*
     * Seek and/or advance to position in the file as described by partition_info
     */
    bool seek_and_fetch_position(File2TileDBBinaryColumnPartitionBase& partition_info, bool& is_read_buffer_exhausted,
        bool force_seek, bool advance_reader);
    /*
     * Return #callsets in current line
     */
    uint64_t get_num_callsets_in_record(const File2TileDBBinaryColumnPartitionBase& partition_info) const
    { return 1u; }
    /*
     * Invoke libcsv function
     * Returns true if (store_in_buffer && buffer is full)
     */
    bool parse_line(const char* line, CSV2TileDBBinaryColumnPartition& csv_partition_info, const unsigned max_token_idx, const bool store_in_buffer);
    /*
     * CSV handling functions - token and end of line
     */
    void handle_token(CSVLineParseStruct* csv_line_parse_ptr, const char* field_ptr, const size_t field_size);
    template<class FieldType>
    void handle_field_token(const char* token_ptr,
        CSVLineParseStruct* csv_line_parse_ptr, CSV2TileDBBinaryColumnPartition& csv_partition_info,
        std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit,
        VariantFieldTypeEnum variant_field_type_enum);
    void handle_end_of_line(CSVLineParseStruct* csv_line_parse_ptr);
  private:
    bool m_cleanup_file;
};

/*
 * Class that maintains state while parsing a csv line
 * Passed to callback function
 */
class CSVLineParseStruct
{
  public:
    CSVLineParseStruct(CSV2TileDBBinary* ptr, CSV2TileDBBinaryColumnPartition* col_info, unsigned max_token_idx,
        bool store_in_buffer)
    {
      m_converter = ptr;
      m_column_partition = col_info;
      m_max_token_idx = max_token_idx;
      m_store_in_buffer = store_in_buffer;
      reset();
    }
    void reset()
    {
      m_token_idx = 0u;
      m_field_idx = 0u;
      m_num_elements_in_field = UINT32_MAX;
      m_field_element_idx = 0u;
      m_row_idx = -1;
      m_enabled_idx_in_file = -1;
      m_cell_size_offset = -1;
    }
    /*
     * Token idx tracking
     */
    void increment_token_idx() { ++m_token_idx; }
    unsigned get_token_idx() const { return m_token_idx; }
    bool is_past_max_token_idx() const { return m_token_idx > m_max_token_idx; }
    /*
     * Field idx tracking
     */
    void increment_field_idx() { ++m_field_idx; }
    unsigned get_field_idx() const { return m_field_idx; }
    /*
     * For fields which are lists of numbers, tracks element in list
     */
    void reset_field_element_idx()
    {
      m_num_elements_in_field = UINT32_MAX;
      m_field_element_idx = 0u;
    }
    bool read_num_elements() const { return (m_num_elements_in_field != UINT32_MAX); }
    void set_num_elements(unsigned val) { m_num_elements_in_field = val; }
    unsigned get_num_elements() const { return m_num_elements_in_field; }
    void increment_field_element_idx()  { ++m_field_element_idx; }
    unsigned get_field_element_idx() const { return m_field_element_idx; }
    /*
     * Access members
     */
    CSV2TileDBBinary* get_csv2tiledb_binary_ptr() { return m_converter; }
    CSV2TileDBBinaryColumnPartition* get_csv_column_partition_ptr() { return m_column_partition; }
    bool should_store_in_buffer() const { return m_store_in_buffer; }
    void set_row_idx(int64_t val) { m_row_idx = val; }
    void set_enabled_idx_in_file(int64_t val) { m_enabled_idx_in_file = val; }
    int64_t get_row_idx() const { return m_row_idx; }
    int64_t get_enabled_idx_in_file() const { return m_enabled_idx_in_file; }
    void set_cell_size_offset(int64_t val) { m_cell_size_offset = val; }
    int64_t get_cell_size_offset() const { return m_cell_size_offset; }
  private:
    bool m_store_in_buffer;
    unsigned m_max_token_idx;
    unsigned m_token_idx;
    unsigned m_field_idx;
    unsigned m_num_elements_in_field;
    unsigned m_field_element_idx;
    CSV2TileDBBinary* m_converter;
    CSV2TileDBBinaryColumnPartition* m_column_partition;
    int64_t m_row_idx;
    int64_t m_enabled_idx_in_file;
    int64_t m_cell_size_offset;
};

#endif
