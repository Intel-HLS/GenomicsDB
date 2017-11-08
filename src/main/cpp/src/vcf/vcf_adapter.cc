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

#ifdef HTSDIR

#include "vcf_adapter.h"
#include "vid_mapper.h"
#include "htslib/tbx.h"

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

//VCFAdapter functions
bool VCFAdapter::add_field_to_hdr_if_missing(bcf_hdr_t* hdr, const VidMapper* id_mapper, const std::string& field_name, int field_type_idx)
{
  auto field_idx = bcf_hdr_id2int(hdr, BCF_DT_ID, field_name.c_str());
  //bcf_hdr_idinfo_exists handles negative field idx
  auto idinfo_exists = bcf_hdr_idinfo_exists(hdr, field_type_idx, field_idx);
  //Field not found
  if(idinfo_exists == 0)
  {
    std::string header_line = "##";
    switch(field_type_idx)
    {
      case BCF_HL_FLT:
        header_line += "FILTER";
        break;
      case BCF_HL_INFO:
        header_line += "INFO";
        break;
      case BCF_HL_FMT:
        header_line += "FORMAT";
        break;
      default:
        throw VCFAdapterException(std::string("Unknown field type ")+std::to_string(field_type_idx));
        break;
    }
    header_line += "=<ID="+field_name;
    if(field_type_idx != BCF_HL_FLT)
    {
      //GT is weird
      if(field_type_idx == BCF_HL_FMT && field_name == "GT")
        header_line += ",Number=1,Type=String,Description=\"Genotype\"";
      else
      {
        assert(id_mapper->get_field_info(field_name));
        auto field_info = *(id_mapper->get_field_info(field_name)); 
        if(field_info.m_bcf_ht_type != BCF_HT_FLAG)
        {
          header_line += ",Number=";
          auto& length_descriptor = field_info.m_length_descriptor;
          if(length_descriptor.get_num_dimensions() == 1u)
          {
            switch(length_descriptor.get_length_descriptor(0u))
            {
              case BCF_VL_FIXED:
                header_line += std::to_string(length_descriptor.get_num_elements());
                break;
              case BCF_VL_VAR:
                header_line += ".";
                break;
              case BCF_VL_A:
                header_line += "A";
                break;
              case BCF_VL_R:
                header_line += "R";
                break;
              case BCF_VL_G:
                header_line += "G";
                break;
              default:
                throw VCFAdapterException("Unhandled field length descriptor "
                    +std::to_string(length_descriptor.get_length_descriptor(0u)));
                break;
            }
          }
        }
        header_line += ",Type=";
        switch(field_info.m_bcf_ht_type)
        {
          case BCF_HT_FLAG:
            header_line += "Flag";
            break;
          case BCF_HT_INT:
            header_line += "Integer";
            break;
          case BCF_HT_REAL:
            header_line += "Float";
            break;
          case BCF_HT_CHAR:
          case BCF_HT_STR:
            header_line += "String";
            break;
          default:
            throw VCFAdapterException("Field type "+std::to_string(field_info.m_bcf_ht_type)+" not handled");
            break;
        }
      }
    }
    header_line += ",Description=\""+field_name+"\"";
    header_line += ">";
    int line_length = 0;
    auto hrec = bcf_hdr_parse_line(hdr, header_line.c_str(), &line_length);
    bcf_hdr_add_hrec(hdr, hrec);
    bcf_hdr_sync(hdr);
    return true;
  }
  else
  {
    const auto* field_info_ptr = id_mapper->get_field_info(field_name);
    assert(field_info_ptr);
    auto field_ht_type = bcf_hdr_id2type(hdr, field_type_idx, field_idx);
    //Don't bother doing any checks for the GT field
    if(field_name != "GT" && field_ht_type != BCF_HT_STR && field_type_idx != BCF_HL_FLT)
    {
      //Allowed configurations - both the JSON and the header specify that:
      //The field is fixed length and agree on the length OR
      //The field is variable length OR
      //field type is BCF_HT_FLAG and VCF header says length is 0 (VCF spec), vid JSON says that length is 1
      if(!((field_info_ptr->m_length_descriptor.is_fixed_length_field()
              && bcf_hdr_id2length(hdr, field_type_idx, field_idx) == BCF_VL_FIXED
              && field_info_ptr->m_length_descriptor.get_num_elements()
              == static_cast<size_t>(bcf_hdr_id2number(hdr, field_type_idx, field_idx)))
            ||
            (!field_info_ptr->m_length_descriptor.is_fixed_length_field()
             && bcf_hdr_id2length(hdr, field_type_idx, field_idx) != BCF_VL_FIXED
            )
            ||
            (field_ht_type == BCF_HT_FLAG
             && field_info_ptr->m_length_descriptor.is_fixed_length_field()
             && field_info_ptr->m_length_descriptor.get_num_elements() == 1
             && bcf_hdr_id2length(hdr, field_type_idx, field_idx) == BCF_VL_FIXED
             && static_cast<int>(bcf_hdr_id2number(hdr, field_type_idx, field_idx)) == 0
            )
          )
        )
        throw VCFAdapterException(std::string("Conflicting field length descriptors and/or field lengths in the vid JSON and VCF header for field ")+field_name);
      //Check for compatible field types
      auto compatible_types = std::unordered_map<int, std::unordered_set<int>>{
           { BCF_HT_FLAG, { BCF_HT_FLAG } },
           { BCF_HT_INT, { BCF_HT_INT } },
           { BCF_HT_REAL, { BCF_HT_REAL } },
           { BCF_HT_INT64, { BCF_HT_INT64 } },
           { BCF_HT_VOID, { BCF_HT_VOID } },
           { BCF_HT_CHAR, { BCF_HT_CHAR, BCF_HT_STR } },
           { BCF_HT_STR, { BCF_HT_CHAR, BCF_HT_STR } }
      };
      if(compatible_types.find(field_ht_type) != compatible_types.end()
          && compatible_types[field_ht_type].find(field_info_ptr->m_bcf_ht_type) == compatible_types[field_ht_type].end())
        throw VCFAdapterException(std::string("Conflicting data types in the vid JSON and VCF header for field ")+field_name);
    }
  }
  return false;
}

VCFAdapter::VCFAdapter(bool open_output, const size_t combined_vcf_records_buffer_size_limit)
{
  m_open_output = open_output;
  m_combined_vcf_records_buffer_size_limit = combined_vcf_records_buffer_size_limit;
  clear();
  m_vcf_header_filename = "";
  m_template_vcf_hdr = 0;
  m_output_fptr = 0;
  m_is_bcf = true;
  m_produce_GT_field = false;
  m_produce_FILTER_field = false;
}

VCFAdapter::~VCFAdapter()
{
  clear();
  if(m_template_vcf_hdr)
    bcf_hdr_destroy(m_template_vcf_hdr);
  if(m_open_output && m_output_fptr)
  {
    bcf_close(m_output_fptr);
    auto status = 0;
    switch(m_index_output_VCF)
    {
      case VCFIndexType::VCF_INDEX_CSI:
        status = bcf_index_build(m_output_filename.c_str(), 14); //bcftools had default 14
        break;
      case VCFIndexType::VCF_INDEX_TBI:
        status = tbx_index_build(m_output_filename.c_str(), 0, &tbx_conf_vcf);
        break;
      default:
        break; //do nothing
    }
    if(status != 0)
      std::cerr << "WARNING: error in creating index for output file "<<m_output_filename<<"\n";
  }
  m_output_fptr = 0;
#ifdef DO_PROFILING
  m_vcf_serialization_timer.print("bcf_t serialization", std::cerr);
#endif
}

void VCFAdapter::clear()
{
  m_reference_genome_info.clear();
  m_vcf_header_filename.clear();
}

void VCFAdapter::initialize(const std::string& reference_genome,
    const std::string& vcf_header_filename,
    std::string output_filename, std::string output_format,
    const size_t combined_vcf_records_buffer_size_limit,
    const bool produce_GT_field,
    const bool index_output_VCF,
    const bool produce_FILTER_field)
{
  //Read template header with fields and contigs
  m_vcf_header_filename = vcf_header_filename;
  if(m_vcf_header_filename.length() > 0u)
  {
    auto* fptr = bcf_open(vcf_header_filename.c_str(), "r");
    m_template_vcf_hdr = bcf_hdr_read_required_sample_line(fptr, 0); //sample line is not required
    bcf_close(fptr);
  }
  else
    m_template_vcf_hdr = initialize_default_header();
  //Output fptr
  std::unordered_map<std::string, bool> valid_output_formats = { {"b", true}, {"bu",true}, {"z",false}, {"",false} };
  if(valid_output_formats.find(output_format) == valid_output_formats.end())
  {
    std::cerr << "INFO: Invalid BCF/VCF output format: "<<output_format<<", will output compressed VCF\n";
    output_format = "z";
  }
  m_is_bcf = valid_output_formats[output_format];
  m_output_filename = output_filename;
  m_index_output_VCF = VCF_INDEX_NONE;
  if(m_open_output)
  {
    m_output_fptr = bcf_open(output_filename.c_str(), ("w"+output_format).c_str());
    if(m_output_fptr == 0)
    {
      std::cerr << "Cannot write to output file "<< output_filename << ", exiting\n";
      exit(-1);
    }
    if(index_output_VCF && !output_filename.empty() && !(output_filename.length() == 1u && output_filename[0] == '-'))
    {
      if(output_format == "z")
        m_index_output_VCF = VCFIndexType::VCF_INDEX_TBI;
      else
      {
        if(output_format == "b")
          m_index_output_VCF = VCFIndexType::VCF_INDEX_CSI;
      }
    }
  }
  //Reference genome
  m_reference_genome_info.initialize(reference_genome);
  m_combined_vcf_records_buffer_size_limit = combined_vcf_records_buffer_size_limit;
  m_produce_GT_field = produce_GT_field;
  m_produce_FILTER_field = produce_FILTER_field;
}

bcf_hdr_t* VCFAdapter::initialize_default_header()
{
  auto hdr = bcf_hdr_init("w");
  bcf_hdr_append(hdr, "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">");
  bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">");
  bcf_hdr_sync(hdr);
  return hdr;
}

void VCFAdapter::print_header()
{
  bcf_hdr_write(m_output_fptr, m_template_vcf_hdr);
}

void VCFAdapter::handoff_output_bcf_line(bcf1_t*& line, const size_t bcf_record_size)
{
  auto write_status = bcf_write(m_output_fptr, m_template_vcf_hdr, line);
  if(write_status != 0)
    throw VCFAdapterException(std::string("Failed to write VCF/BCF record at position ")
        +bcf_hdr_id2name(m_template_vcf_hdr, line->rid)+", "
        +std::to_string(line->pos+1));
}

BufferedVCFAdapter::BufferedVCFAdapter(unsigned num_circular_buffers, unsigned max_num_entries, const size_t combined_vcf_records_buffer_size_limit)
  : VCFAdapter(true, combined_vcf_records_buffer_size_limit), CircularBufferController(num_circular_buffers)
{
  clear();
  m_line_buffers.resize(num_circular_buffers);
  m_num_valid_entries.resize(num_circular_buffers, 0u);
  m_combined_vcf_records_buffer_sizes.resize(num_circular_buffers, 0ull);
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
  m_combined_vcf_records_buffer_sizes.clear();
}

void BufferedVCFAdapter::handoff_output_bcf_line(bcf1_t*& line, const size_t bcf_record_size)
{
  auto write_idx = get_write_idx();
  auto& line_buffer = m_line_buffers[write_idx];
  //Need to resize buffer - non-common case
  if(m_num_valid_entries[write_idx] >= line_buffer.size())
    resize_line_buffer(line_buffer, 2u*line_buffer.size()+1u);
  assert(m_num_valid_entries[write_idx] < line_buffer.size());
  std::swap<bcf1_t*>(line, line_buffer[m_num_valid_entries[write_idx]]);
  ++(m_num_valid_entries[write_idx]);
  m_combined_vcf_records_buffer_sizes[write_idx] += bcf_record_size;
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
  auto write_idx = get_write_idx();
  //Advance if something was written
  if(m_num_valid_entries[write_idx])
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
    auto write_status = bcf_write(m_output_fptr, m_template_vcf_hdr, m_line_buffers[read_idx][i]);
    if(write_status != 0)
      throw VCFAdapterException(std::string("Failed to write VCF/BCF record at position ")
          +bcf_hdr_id2name(m_template_vcf_hdr, m_line_buffers[read_idx][i]->rid)+", "
          +std::to_string(m_line_buffers[read_idx][i]->pos+1));
  }
  m_num_valid_entries[read_idx] = 0u;
  m_combined_vcf_records_buffer_sizes[read_idx] = 0ull;
  advance_read_idx();
}

void VCFSerializedBufferAdapter::print_header()
{
  assert(m_rw_buffer);
  auto offset = bcf_hdr_serialize(m_template_vcf_hdr, &(m_rw_buffer->m_buffer[0]), m_rw_buffer->m_num_valid_bytes, m_rw_buffer->m_buffer.size(),
      m_is_bcf ? 1u : 0u, m_keep_idx_fields_in_bcf_header ? 1u : 0u);
  //Buffer capacity was too small, resize
  while(offset == m_rw_buffer->m_num_valid_bytes)
  {
    m_rw_buffer->m_buffer.resize(2u*(m_rw_buffer->m_buffer.size())+1u);
    offset = bcf_hdr_serialize(m_template_vcf_hdr, &(m_rw_buffer->m_buffer[0]), m_rw_buffer->m_num_valid_bytes, m_rw_buffer->m_buffer.size(),
        m_is_bcf ? 1u : 0u, m_keep_idx_fields_in_bcf_header ? 1u : 0u);
  }
  m_rw_buffer->m_num_valid_bytes = offset;
}

void VCFSerializedBufferAdapter::handoff_output_bcf_line(bcf1_t*& line, const size_t bcf_record_size)
{
#ifdef DO_PROFILING
  m_vcf_serialization_timer.start();
#endif
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
#ifdef DO_PROFILING
  m_vcf_serialization_timer.stop();
#endif
}

#endif //ifdef HTSDIR
