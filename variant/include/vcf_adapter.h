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
    ~VCFAdapter();
    void clear();
    void initialize(const std::string& reference_genome, const std::string& vcf_header_filename,
        std::string output_filename, std::string output_format="");
    bcf_hdr_t* get_vcf_header() { return m_template_vcf_hdr; }
    void print_bcf_line(bcf1_t* line) { bcf_write(m_output_fptr, m_template_vcf_hdr, line); }
    void print_header();
    char get_reference_base_at_position(const char* contig, int pos)
    { return m_reference_genome_info.get_reference_base_at_position(contig, pos); }
  private:
    //Template VCF header to start with
    std::string m_vcf_header_filename;
    bcf_hdr_t* m_template_vcf_hdr;
    //Reference genome info
    ReferenceGenomeInfo m_reference_genome_info;
    //Output fptr
    htsFile* m_output_fptr;
    bool m_is_bcf;
};

#endif  //ifdef HTSDIR

#endif
