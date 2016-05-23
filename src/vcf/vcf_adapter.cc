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

#ifdef HTSDIR

#include "vcf_adapter.h"
#include "vid_mapper.h"

//ReferenceGenomeInfo functions
void ReferenceGenomeInfo::initialize(const std::string& reference_genome)
{
  m_reference_faidx = fai_load(reference_genome.c_str());
  assert(m_reference_faidx);
  m_reference_last_seq_read = "";
  //buffer
  m_buffer.resize(32768u+8);     //32KB
}

char ReferenceGenomeInfo::get_reference_base_at_position(const char* contig, int pos)
{
  //See if pos is within the last buffer read
  if(strcmp(m_reference_last_seq_read.c_str(), contig) == 0 && m_reference_last_read_pos <= pos)
  {
    int offset = pos - m_reference_last_read_pos;
    if(offset < m_reference_num_bases_read)
      return m_buffer[offset];
  }
  int length = 0;
  faidx_fetch_seq_into_buffer(m_reference_faidx, contig, pos, pos+m_buffer.size()-8u, &(m_buffer[0]), &length);
  assert(length > 0);
  m_reference_last_seq_read = contig;
  m_reference_last_read_pos = pos;
  m_reference_num_bases_read = length;
  return m_buffer[0];
}

VCFAdapter::VCFAdapter()
{
  clear();
  m_vcf_header_filename = "";
  m_template_vcf_hdr = 0;
  m_output_fptr = 0;
  m_is_bcf = true;
}

VCFAdapter::~VCFAdapter()
{
  clear();
  if(m_template_vcf_hdr)
    bcf_hdr_destroy(m_template_vcf_hdr);
  if(m_output_fptr)
    bcf_close(m_output_fptr);
}

void VCFAdapter::clear()
{
  m_reference_genome_info.clear();
  m_vcf_header_filename.clear();
}

void VCFAdapter::initialize(const std::string& reference_genome,
    const std::string& vcf_header_filename,
    std::string output_filename, std::string output_format)
{
  //Read template header with fields and contigs
  m_vcf_header_filename = vcf_header_filename;
  auto* fptr = bcf_open(vcf_header_filename.c_str(), "r");
  m_template_vcf_hdr = bcf_hdr_read(fptr);
  bcf_close(fptr);
  //Output fptr
  std::unordered_map<std::string, bool> valid_output_formats = { {"b", true}, {"bu",true}, {"z",false}, {"",false} };
  if(valid_output_formats.find(output_format) == valid_output_formats.end())
  {
    std::cerr << "INFO: Invalid BCF/VCF output format: "<<output_format<<", will output compressed VCF\n";
    output_format = "z";
  }
  m_is_bcf = valid_output_formats[output_format];
  m_output_fptr = bcf_open(output_filename.c_str(), ("w"+output_format).c_str());
  m_output_filename = output_filename;
  if(m_output_fptr == 0)
  {
    std::cerr << "Cannot write to output file "<< output_filename << ", exiting\n";
    exit(-1);
  }
  //Reference genome
  m_reference_genome_info.initialize(reference_genome);
}

void VCFAdapter::print_header()
{
  bcf_hdr_write(m_output_fptr, m_template_vcf_hdr);
}

BufferedVCFAdapter::BufferedVCFAdapter(unsigned num_circular_buffers, unsigned max_num_entries)
  : VCFAdapter(), CircularBufferController(num_circular_buffers)
{
  clear();
  m_line_buffers.resize(num_circular_buffers);
  m_num_valid_entries.resize(num_circular_buffers);
  for(auto i=0u;i<m_num_valid_entries.size();++i)
    m_num_valid_entries[i] = 0u;
  //Initialize buffers
  for(auto& line_buffer : m_line_buffers)
    resize_line_buffer(line_buffer, max_num_entries);
}

BufferedVCFAdapter::~BufferedVCFAdapter()
{
  for(auto& line_buffer : m_line_buffers)
    for(auto& line : line_buffer)
      bcf_destroy(line);
  clear();
}

void BufferedVCFAdapter::clear()
{
  m_line_buffers.clear();
  m_num_valid_entries.clear();
}

void BufferedVCFAdapter::handoff_output_bcf_line(bcf1_t*& line)
{
  auto write_idx = get_write_idx();
  auto& line_buffer = m_line_buffers[write_idx];
  //Need to resize buffer - non-common case
  if(m_num_valid_entries[write_idx] >= line_buffer.size())
    resize_line_buffer(line_buffer, 2u*line_buffer.size()+1u);
  assert(m_num_valid_entries[write_idx] < line_buffer.size());
  std::swap<bcf1_t*>(line, line_buffer[m_num_valid_entries[write_idx]]);
  ++(m_num_valid_entries[write_idx]);
}

void BufferedVCFAdapter::resize_line_buffer(std::vector<bcf1_t*>& line_buffer, unsigned new_size)
{
  if(new_size <= line_buffer.size())     //never reduce
    return;
  auto curr_idx = line_buffer.size();
  line_buffer.resize(new_size);
  for(auto i=curr_idx;i<line_buffer.size();++i)
    line_buffer[i] = bcf_init();
}

void BufferedVCFAdapter::advance_write_idx()
{
  //Advance if something was written
  if(m_num_valid_entries[get_write_idx()])
    CircularBufferController::advance_write_idx();
}

void BufferedVCFAdapter::do_output()
{
  if(get_num_entries_with_valid_data() == 0u)
    return;
  auto read_idx = get_read_idx();
  assert(m_num_valid_entries[read_idx] <= m_line_buffers[read_idx].size());
  for(auto i=0u;i<m_num_valid_entries[read_idx];++i)
  {
    assert(m_line_buffers[read_idx][i]);
    bcf_write(m_output_fptr, m_template_vcf_hdr, m_line_buffers[read_idx][i]);
  }
  m_num_valid_entries[read_idx] = 0u;
  advance_read_idx();
}

void VCFSerializedBufferAdapter::print_header()
{
  assert(m_rw_buffer);
  auto offset = bcf_hdr_serialize(m_template_vcf_hdr, &(m_rw_buffer->m_buffer[0]), m_rw_buffer->m_num_valid_bytes, m_rw_buffer->m_buffer.size(),
      m_is_bcf ? 1u : 0u);
  //Buffer capacity was too small, resize
  while(offset == m_rw_buffer->m_num_valid_bytes)
  {
    m_rw_buffer->m_buffer.resize(2u*(m_rw_buffer->m_buffer.size())+1u);
    offset = bcf_hdr_serialize(m_template_vcf_hdr, &(m_rw_buffer->m_buffer[0]), m_rw_buffer->m_num_valid_bytes, m_rw_buffer->m_buffer.size(),
       m_is_bcf ? 1u : 0u);
  }
  m_rw_buffer->m_num_valid_bytes = offset;
}

void VCFSerializedBufferAdapter::handoff_output_bcf_line(bcf1_t*& line)
{
  assert(m_rw_buffer);
  auto offset = bcf_serialize(line, &(m_rw_buffer->m_buffer[0]), m_rw_buffer->m_num_valid_bytes, m_rw_buffer->m_buffer.size(),
       m_is_bcf ? 1u : 0u, m_template_vcf_hdr, &m_hts_string);
  //Buffer capacity was too small, resize
  while(offset == m_rw_buffer->m_num_valid_bytes)
  {
    m_rw_buffer->m_buffer.resize(2u*(m_rw_buffer->m_buffer.size())+1u);
    offset = bcf_serialize(line, &(m_rw_buffer->m_buffer[0]), m_rw_buffer->m_num_valid_bytes, m_rw_buffer->m_buffer.size(),
       m_is_bcf ? 1u : 0u, m_template_vcf_hdr, &m_hts_string);
  }
  m_rw_buffer->m_num_valid_bytes = offset;
}

#endif //ifdef HTSDIR
