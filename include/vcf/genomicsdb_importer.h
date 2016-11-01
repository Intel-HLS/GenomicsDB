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
#ifndef GENOMICSDB_IMPORTER_H
#define GENOMICSDB_IMPORTER_H

#include "tiledb_loader.h"

class GenomicsDBImporterException : public std::exception {
  public:
    GenomicsDBImporterException(const std::string m="") : msg_("GenomicsDBImporterException : "+m) { ; }
    ~GenomicsDBImporterException() { ; }
    // ACCESSORS
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class GenomicsDBImporter
{
  public:
    GenomicsDBImporter(const std::string& loader_config_file, const int rank,
        const int64_t lb_callset_row_idx=0, const int64_t ub_callset_row_idx=INT64_MAX-1)
    {
      m_is_loader_setup = false;
      m_rank = rank;
      m_loader_config_file = loader_config_file;
      m_lb_callset_row_idx = lb_callset_row_idx;
      m_ub_callset_row_idx = ub_callset_row_idx;
      m_loader_ptr = 0;
      m_read_state = 0;
    }
    //Delete copy constructor
    GenomicsDBImporter(const GenomicsDBImporter& other) = delete;
    //Define move constructor
    GenomicsDBImporter(GenomicsDBImporter&& other);
    //Destructor
    ~GenomicsDBImporter();
    //Buffer streams
    void add_buffer_stream(const std::string& name, const VidFileTypeEnum buffer_stream_type,
        const size_t capacity)
    {
      add_buffer_stream(name, buffer_stream_type, capacity, 0, 0);
    }
    void add_buffer_stream(const std::string& name, const VidFileTypeEnum buffer_stream_type,
        const size_t capacity, const uint8_t* initialization_buffer, const size_t num_bytes_in_initialization_buffer);
    /*
     * No more buffer streams can be added to this object after setup_loader is called
     */
    void setup_loader();
    /*
     * Can be used by callers to pre-allocate vector<BufferStreamIdentifier>
     */
    size_t get_max_num_buffer_stream_identifiers() const { return m_loader_ptr->get_max_num_buffer_stream_identifiers(); }
    /*
     * Write to buffer stream
     */
    void write_data_to_buffer_stream(const int64_t buffer_stream_idx, const unsigned partition_idx, const uint8_t* data, const size_t num_bytes)
    {
      if(!m_is_loader_setup)
        throw GenomicsDBImporterException("Cannot write data to buffer stream in the GenomicsDBImporter without calling setup_loader() first");
      assert(m_loader_ptr);
      m_loader_ptr->write_data_to_buffer_stream(buffer_stream_idx, partition_idx, data, num_bytes);
    }
    /*
     * Import next batch of data into TileDB/GenomicsDB
     */
    void import_batch();
    /*
     * Obtains buffer stream identifiers that are exhausted and must be replenished by the caller
     */
    const std::vector<BufferStreamIdentifier>& get_exhausted_buffer_stream_identifiers() const
    {
      return m_loader_ptr->get_exhausted_buffer_stream_identifiers();
    }
    bool is_done() const { return m_read_state->is_done(); }
    void finish()
    {
      m_loader_ptr->finish_read_all(*m_read_state);
    }
  private:
    void copy_simple_members(const GenomicsDBImporter& other);
  private:
    bool m_is_loader_setup;
    int m_rank;
    std::string m_loader_config_file;
    int64_t m_lb_callset_row_idx;
    int64_t m_ub_callset_row_idx;
    std::vector<BufferStreamInfo> m_buffer_stream_info_vec;
    VCF2TileDBLoader* m_loader_ptr;
    VCF2TileDBLoaderReadState* m_read_state;
};

#endif
