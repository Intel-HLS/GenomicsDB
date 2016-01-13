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

#endif //ifdef HTSDIR
