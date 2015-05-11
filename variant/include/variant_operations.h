#ifndef VARIANT_OPERATIONS_H
#define VARIANT_OPERATIONS_H

#include "variant.h"
#include "lut.h"

typedef void (*scan_operator_type)(Variant& variant, void* data);

//Base wrapper class for any data that requires re-mapping based on allele order
//Will be used as argument to remap_ functions
//Tthe goal of this class (actually its subclasses) is to return an address 
//at which the remap_ functions will store data for a given call_idx and the allele/genotype
//index. The remap_ functions are responsible for casting the pointer to the correct type
class RemappedDataWrapperBase
{
  public:
    virtual void* put_address(uint64_t input_call_idx, unsigned allele_or_gt_idx) = 0; 
};

/*
 * Subclass that redirects output of the remap_ functions to a matrix of values. The
 * pointer returned by the put function is an address into the matrix
 * The matrix is organized in "genotype major" order, ie, the rows correspond to 
 * genotypes while columns correspond to input calls
 */
template<class DataType>
class RemappedMatrix : public RemappedDataWrapperBase
{
  public:
    virtual void* put_address(uint64_t input_call_idx, unsigned allele_or_gt_idx);
    void resize(uint64_t num_rows, uint64_t num_columns, DataType init_value);
    std::vector<std::vector<DataType>>& get() { return m_matrix; }
    const std::vector<std::vector<DataType>>& get() const { return m_matrix; }
  private:
    std::vector<std::vector<DataType>> m_matrix;
};

/*
 * Class that helps remap PL and other allele dependent fields within the Variant object itself
 * The address returned is the address within the VariantField object corresponding to element allele_or_gt_idx
 * for the call corresponding to input_call_idx.
 * The VariantFieldObject is determined by the m_queried_field_idx value
 */
class RemappedVariant : public RemappedDataWrapperBase
{
  public:
    RemappedVariant(Variant& variant, unsigned queried_field_idx)
      : RemappedDataWrapperBase(), m_variant(&variant), m_queried_field_idx(queried_field_idx)
    {}
    virtual void* put_address(uint64_t input_call_idx, unsigned allele_or_gt_idx);
  private:
    unsigned m_queried_field_idx;
    Variant* m_variant;
};

class VariantOperations
{
  public:
    /*
     * Obtains a merged REF field as defined in BCF spec
     */
    static void merge_reference_allele(const Variant& variant, const VariantQueryConfig& query_config,
        std::string& merged_reference_allele);
    /*
     * Obtains a merged ALT list as defined in BCF spec
     */
    static const std::vector<std::string>  merge_alt_alleles(const Variant& variant,
        const VariantQueryConfig& query_config,
        const std::string& merged_reference_allele,
        CombineAllelesLUT& alleles_LUT, bool& NON_REF_exists);
    /*
     * Reorders fields whose length and order depend on the number of alleles (BCF_VL_R or BCF_VL_A)
     */
    template<class DataType>
    static void remap_data_based_on_alleles(const std::vector<DataType>& input_data,
        const uint64_t input_call_idx, 
        const CombineAllelesLUT& alleles_LUT, const unsigned num_merged_alleles, bool NON_REF_exists, bool alt_alleles_only,
        RemappedDataWrapperBase& remapped_data,
        std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value);
    /*
     * Reorders fields whose length and order depend on the number of genotypes (BCF_VL_G)
     */
    template<class DataType>
    static void remap_data_based_on_genotype(const std::vector<DataType>& input_data,
        const uint64_t input_call_idx, 
        const CombineAllelesLUT& alleles_LUT, const unsigned num_merged_alleles, bool NON_REF_exists,
        RemappedDataWrapperBase& remapped_data,
        std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value);
    static void do_dummy_genotyping(Variant& variant, std::ostream& output);
};

/*
 * Base class for all operations that want to operate on Variant objects created by scan/gt_get_column functions
 */
class SingleVariantOperatorBase
{
  public:
    SingleVariantOperatorBase()
    {
      clear();
      m_NON_REF_exists = false;
    }
    void clear();
    /*
     * Main operate function - should be overridden by all child classes
     * Basic operate creates:
     * (a) Merged reference allele
     * (b) Merged ALT allele list and updates the alleles LUT
     */
    virtual void operate(Variant& variant, const VariantQueryConfig& query_config);
  protected:
    //Maintain mapping between alleles in input VariantCalls and merged allele list
    CombineAllelesLUT m_alleles_LUT;
    /*Flag that determines whether the merged ALT list contains NON_REF allele*/
    bool m_NON_REF_exists;
    //Merged reference allele and ALT alleles
    std::string m_merged_reference_allele;
    std::vector<std::string> m_merged_alt_alleles;
};

class DummyGenotypingOperator : public SingleVariantOperatorBase
{
  public:
    DummyGenotypingOperator() : SingleVariantOperatorBase() { m_output_stream = &(std::cout); }
    virtual void operate(Variant& variant, const VariantQueryConfig& query_config);
    std::ostream* m_output_stream;
};

/*
 * Copies info in Variant object into its result vector
 */
class GA4GHOperator : public SingleVariantOperatorBase
{
  public:
    void clear() { m_variants.clear(); }
    virtual void operate(Variant& variant, const VariantQueryConfig& query_config);
    const std::vector<Variant>& get_variants() const { return m_variants; }
    std::vector<Variant>& get_variants() { return m_variants; }
  private:
    std::vector<Variant> m_variants;
};

#endif
