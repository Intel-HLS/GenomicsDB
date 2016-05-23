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

#ifndef VCF_ADAPTER_H
#define VCF_ADAPTER_H

#ifdef HTSDIR

#include "headers.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"

//Struct for accessing reference genome
//reference genome access - required for gVCF merge
class ReferenceGenomeInfo
{
  public:
    ReferenceGenomeInfo()
    {
      clear();
      m_reference_faidx = 0;
      m_reference_last_read_pos = -1;
      m_reference_num_bases_read = 0;
      m_reference_last_seq_read = "";
    }
    void clear()
    {
      m_reference_last_seq_read.clear();
      m_buffer.clear();
    }
    ~ReferenceGenomeInfo()
    {
      clear();
      if(m_reference_faidx)
        fai_destroy(m_reference_faidx);
    }
    void initialize(const std::string& reference_genome);
    char get_reference_base_at_position(const char* contig, int pos);
  private:
    int m_reference_last_read_pos;
    int m_reference_num_bases_read;
    std::string m_reference_last_seq_read;
    std::vector<char> m_buffer;
    faidx_t* m_reference_faidx;
};

class VCFAdapter
{
  public:
    VCFAdapter();
    virtual ~VCFAdapter();
    void clear();
    void initialize(const std::string& reference_genome, const std::string& vcf_header_filename,
        std::string output_filename, std::string output_format="");
    bcf_hdr_t* get_vcf_header() { return m_template_vcf_hdr; }
    /*
     * The line is ready for output
     * Child classes might actually just swap out the pointer so that the actual output is performed by
     * a thread off the critical path
     **/
    virtual void handoff_output_bcf_line(bcf1_t*& line) { bcf_write(m_output_fptr, m_template_vcf_hdr, line); }
    virtual void print_header();
    /*
     * Return true in child class if some output causes buffer to be full. Default: return false
     */
    virtual bool overflow() const { return false; }
    char get_reference_base_at_position(const char* contig, int pos)
    { return m_reference_genome_info.get_reference_base_at_position(contig, pos); }
  protected:
    //Output file
    std::string m_output_filename;
    //Template VCF header to start with
    std::string m_vcf_header_filename;
    bcf_hdr_t* m_template_vcf_hdr;
    //Reference genome info
    ReferenceGenomeInfo m_reference_genome_info;
    //Output fptr
    htsFile* m_output_fptr;
    bool m_is_bcf;
};

class BufferedVCFAdapter : public VCFAdapter, public CircularBufferController
{
  public:
    BufferedVCFAdapter(unsigned num_circular_buffers, unsigned max_num_entries);
    virtual ~BufferedVCFAdapter();
    void clear();
    virtual void handoff_output_bcf_line(bcf1_t*& line);
    void advance_write_idx();
    void do_output();
  private:
    void resize_line_buffer(std::vector<bcf1_t*>& line_buffer, unsigned new_size);
    std::vector<std::vector<bcf1_t*>> m_line_buffers;   //Outer vector for double-buffering
    std::vector<unsigned> m_num_valid_entries;  //One per double-buffer
};

class VCFSerializedBufferAdapter: public VCFAdapter
{
  public:
    VCFSerializedBufferAdapter(const size_t overflow_limit, bool print_output)
      : VCFAdapter()
    {
      m_rw_buffer = 0;
      m_overflow_limit = overflow_limit;
      //Temporary hts string
      m_hts_string.l = 0u;
      m_hts_string.m = 4096u;
      m_hts_string.s = (char*)malloc(m_hts_string.m);
      assert(m_hts_string.s);
      m_write_fptr = print_output ? (m_output_filename == "" ? stdout : fopen(m_output_filename.c_str(), "w")) : 0;
    }
    ~VCFSerializedBufferAdapter()
    {
      if(m_write_fptr && m_write_fptr != stdout)
        fclose(m_write_fptr);
      m_write_fptr = 0;
    }
    void set_buffer(RWBuffer& buffer) { m_rw_buffer = &buffer; }
    void print_header();
    void handoff_output_bcf_line(bcf1_t*& line);
    inline bool overflow() const
    {
      assert(m_rw_buffer);
      return (m_rw_buffer->m_num_valid_bytes >= m_overflow_limit);
    }
    void do_output()
    {
      assert(m_write_fptr);
      assert(m_rw_buffer);
      fwrite(&(m_rw_buffer->m_buffer[0]), 1u,  m_rw_buffer->m_num_valid_bytes, m_write_fptr);
    }
  private:
    RWBuffer* m_rw_buffer;
    FILE* m_write_fptr;
    kstring_t m_hts_string;
    size_t m_overflow_limit;
};

#endif  //ifdef HTSDIR

#endif
