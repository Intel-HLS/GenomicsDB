#ifndef VARIANT_OPERATIONS_H
#define VARIANT_OPERATIONS_H

#include "variant.h"
#include "lut.h"

typedef void (*scan_operator_type)(const GTColumn*, void* data);

class VariantOperations
{
  public:
    static void merge_reference_allele(const std::vector<std::string>& reference_vector, std::string& merged_reference_allele);
    static const std::vector<std::string>  merge_alt_alleles(const std::vector<std::string>& reference_vector, 
            const std::vector<std::vector<std::string>> alt_alleles,
            const std::string& merged_reference_allele,
            CombineAllelesLUT& alleles_LUT);
    static void  remap_pl(const std::vector<int32_t>& pls, const uint32_t input_sample_idx,
            const CombineAllelesLUT& alleles_LUT, const unsigned num_merged_alleles,
            std::vector<std::vector<int>>& remapped_pls,
            std::vector<uint64_t>& num_missing_samples);
    static void do_dummy_genotyping(const GTColumn* gt_column, std::ostream& output);
};

#endif
