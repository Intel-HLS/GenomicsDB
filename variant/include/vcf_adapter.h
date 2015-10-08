#ifndef VCF_ADAPTER_H
#define VCF_ADAPTER_H

#include "headers.h"
#include "htslib/vcf.h"
#include "bcftools.h"

class VCFAdapter
{
  public:
    VCFAdapter();
    ~VCFAdapter();
    void clear();
    void initialize(const std::string& sqlite_filename, const std::string& vcf_header_filename);
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
  private:
    //Template VCF header to start with
    std::string m_vcf_header_filename;
    bcf_hdr_t* m_template_vcf_hdr;
    //Mapping from sample names to idx
    std::string m_sqlite_filename;
    sqlite_mappings_struct m_sqlite_mapping_info;
    //sorted vectors of pair<contig begin/end, idx> 
    std::vector<std::pair<int64_t, int>> m_contig_begin_2_idx;
    std::vector<std::pair<int64_t, int>> m_contig_end_2_idx;
};

#endif
