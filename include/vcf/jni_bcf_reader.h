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

#ifndef JNI_BCF_READER_H
#define JNI_BCF_READER_H

#include "headers.h"
#include "broad_combined_gvcf.h"
#include "variant_storage_manager.h"
#include "query_variants.h"

class GenomicsDBJNIException : public std::exception {
  public:
    GenomicsDBJNIException(const std::string m="") : msg_("GenomicsDBJNIException : "+m) { ; } 
    ~GenomicsDBJNIException() { ; } 
    // ACCESSORS
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class JNIBCFReader
{
  public:
    JNIBCFReader(const std::string& loader_config_file, const std::string& query_config_file,
        const char* chr, const int start, const int end,
        int my_rank=0, size_t buffer_capacity=1048576u, size_t tiledb_segment_size=1048576u, const char* output_format="bu");
    JNIBCFReader(const std::string& loader_config_file, const std::string& query_config_file, int my_rank=0,
        size_t buffer_capacity=1048576u, size_t tiledb_segment_size=1048576u, const char* output_format="bu")
      : JNIBCFReader(loader_config_file, query_config_file,
          "", 0, 0,
          my_rank, buffer_capacity, tiledb_segment_size, output_format)
      {  }
    //Delete copy and move constructors
    JNIBCFReader(const JNIBCFReader& other) = delete;
    JNIBCFReader(JNIBCFReader&& other) = delete;
    ~JNIBCFReader();
    size_t get_buffer_capacity() const
    {
      return m_buffers[0u].m_buffer.size();
    }
    const RWBuffer& get_read_batch() const
    {
      return m_buffers[m_buffer_control.get_read_idx()];
    }
    size_t read_and_advance(uint8_t* dst, size_t offset, size_t n);
    void produce_next_batch();
    uint8_t read_next_byte();
    inline bool end() const { return m_done; }
  private:
    void set_write_buffer();
    void reset_read_buffer();
  private:
    bool m_done;
    FileBasedVidMapper m_vid_mapper;
    VariantStorageManager* m_storage_manager;
    VariantQueryProcessor* m_query_processor;
    VariantQueryConfig m_query_config;
    VCFSerializedBufferAdapter m_vcf_adapter;
    unsigned m_query_column_interval_idx;
    VariantQueryProcessorScanState m_scan_state;
    BroadCombinedGVCFOperator* m_combined_bcf_operator;
    //If using ping-pong buffering, then multiple buffers exist
    std::vector<RWBuffer> m_buffers;
    CircularBufferController m_buffer_control;
};

#endif
