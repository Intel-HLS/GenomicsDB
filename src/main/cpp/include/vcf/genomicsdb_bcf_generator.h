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

#ifndef GENOMICSDB_BCF_GENERATOR_H
#define GENOMICSDB_BCF_GENERATOR_H

#include "headers.h"
#include "broad_combined_gvcf.h"
#include "variant_storage_manager.h"
#include "query_variants.h"
#include "timer.h"
#include "genomicsdb_jni_exception.h"

class GenomicsDBBCFGenerator
{
  public:
    GenomicsDBBCFGenerator(const std::string& loader_config_file, const std::string& query_config_file,
        const char* chr, const int start, const int end,
        int my_rank=0, size_t buffer_capacity=DEFAULT_COMBINED_VCF_RECORDS_BUFFER_SIZE, size_t tiledb_segment_size=1048576u,
        const char* output_format="bu",
        const bool produce_header_only=false,
        const bool use_missing_values_only_not_vector_end=false,
        const bool keep_idx_fields_in_bcf_header=true);
    GenomicsDBBCFGenerator(const std::string& loader_config_file, const std::string& query_config_file, int my_rank=0,
        size_t buffer_capacity=DEFAULT_COMBINED_VCF_RECORDS_BUFFER_SIZE, size_t tiledb_segment_size=1048576u, const char* output_format="bu",
        const bool produce_header_only=false,
        const bool use_missing_values_only_not_vector_end=false)
      : GenomicsDBBCFGenerator(loader_config_file, query_config_file,
          "", 0, 0,
          my_rank, buffer_capacity, tiledb_segment_size, output_format,
          produce_header_only,
          use_missing_values_only_not_vector_end)
      {  }
    //Delete copy and move constructors
    GenomicsDBBCFGenerator(const GenomicsDBBCFGenerator& other) = delete;
    GenomicsDBBCFGenerator(GenomicsDBBCFGenerator&& other) = delete;
    ~GenomicsDBBCFGenerator();
    size_t get_buffer_capacity() const
    {
      return m_buffers[0u].m_buffer.size();
    }
    const RWBuffer& get_read_batch() const
    {
      return m_buffers[m_buffer_control.get_read_idx()];
    }
    /*
     * if n == SIZE_MAX, simply produce the next batch, but don't advance pointers within buffer
     * if dst not NULL, copy n bytes to dst+offset
     */
    size_t read_and_advance(uint8_t* dst, size_t offset, size_t n);
    uint8_t read_next_byte();
    inline bool end() const { return m_done; }
  private:
    void set_write_buffer();
    void reset_read_buffer();
    void produce_next_batch();
  private:
    bool m_done;
    bool m_produce_header_only;
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
#ifdef DO_PROFILING
    Timer m_timer;
#endif
};

#endif
