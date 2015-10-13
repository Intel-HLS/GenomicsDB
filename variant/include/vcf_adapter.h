#ifndef VCF_ADAPTER_H
#define VCF_ADAPTER_H

#ifdef HTSDIR

#include "headers.h"
#include "htslib/vcf.h"
#include "bcftools.h"

class VCFAdapter
{
  public:
    VCFAdapter();
    ~VCFAdapter();
    void clear();
    void initialize(const std::string& sqlite_filename, const std::string& reference_genome,
        const std::string& vcf_header_filename,
        std::string output_filename, std::string output_format="");
    /*
     * Given a position in a flattened 'address' space [TileDB column idx], get the contig_name and location
     * in the contig [0-based]
     * Returns true if valid contig found, false otherwise
     */
    bool get_contig_location(int64_t position, std::string& contig_name, int64_t& contig_position) const;
    /*
     * Given a position in a flattened 'address' space [TileDB column idx], get the next contig_name and starting
     * location of the contig in the flattened space 
     * Returns true if valid contig found, false otherwise
     */
    bool get_next_contig_location(int64_t position, std::string& next_contig_name, int64_t& next_contig_offset) const;
    bcf_hdr_t* get_vcf_header() { return m_template_vcf_hdr; }
    void print_bcf_line(bcf1_t* line) { bcf_write(m_output_fptr, m_template_vcf_hdr, line); }
    void print_header();
    const char* get_sample_name_for_idx(int64_t sample_idx) const;
    char get_reference_base_at_position(const char* contig, int pos);
  private:
    //Template VCF header to start with
    std::string m_vcf_header_filename;
    bcf_hdr_t* m_template_vcf_hdr;
    //Reference genome info
    reference_genome_info m_reference_genome_info;
    //Mapping from sample names to idx
    std::string m_sqlite_filename;
    sqlite_mappings_struct m_sqlite_mapping_info;
    //sorted vectors of pair<contig begin/end, idx> 
    std::vector<std::pair<int64_t, int>> m_contig_begin_2_idx;
    std::vector<std::pair<int64_t, int>> m_contig_end_2_idx;
    //Output fptr
    htsFile* m_output_fptr;
    bool m_is_bcf;
};

#endif  //ifdef HTSDIR

#endif
