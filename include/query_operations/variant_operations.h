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
    Variant* m_variant;
    unsigned m_queried_field_idx;
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
    static void merge_alt_alleles(const Variant& variant,
        const VariantQueryConfig& query_config,
        const std::string& merged_reference_allele,
        CombineAllelesLUT& alleles_LUT, std::vector<std::string>& merged_alt_alleles, bool& NON_REF_exists);
    /*
     * Remaps GT field of Calls in the combined Variant based on new allele order
     */
    static void remap_GT_field(const std::vector<int>& input_GT, std::vector<int>& output_GT,
        const CombineAllelesLUT& alleles_LUT, const uint64_t input_call_idx);
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
      m_remapping_needed = true;
    }
    void clear();
    /*
     * Main operate function - should be overridden by all child classes
     * Basic operate creates:
     * (a) Merged reference allele
     * (b) Merged ALT allele list and updates the alleles LUT
     */
    virtual void operate(Variant& variant, const VariantQueryConfig& query_config);
    /*
     * Return true in child class if some output buffer used by the operator
     * is full. Default implementation: return false
     */
    virtual bool overflow() const { return false; }
  protected:
    //Maintain mapping between alleles in input VariantCalls and merged allele list
    CombineAllelesLUT m_alleles_LUT;
    /*Flag that determines whether the merged ALT list contains NON_REF allele*/
    bool m_NON_REF_exists;
    //Merged reference allele and ALT alleles
    std::string m_merged_reference_allele;
    std::vector<std::string> m_merged_alt_alleles;
    //Flag that determines if any allele re-ordering occurred and whether fields such
    //as PL/AD need to be re-ordered
    bool m_remapping_needed;
};

class DummyGenotypingOperator : public SingleVariantOperatorBase
{
  public:
    DummyGenotypingOperator() : SingleVariantOperatorBase() { m_output_stream = &(std::cout); }
    virtual void operate(Variant& variant, const VariantQueryConfig& query_config);
    std::ostream* m_output_stream;
};

template<class T>
T get_zero_value() { return 0; }

//Base virtual class for storing a big bag of handler functions
class VariantFieldHandlerBase 
{
  public:
    VariantFieldHandlerBase() { ; }
    virtual ~VariantFieldHandlerBase() = default;
    virtual void remap_vector_data(std::unique_ptr<VariantFieldBase>& orig_field_ptr, uint64_t curr_call_idx_in_variant, 
        const CombineAllelesLUT& alleles_LUT, unsigned num_merged_alleles, bool non_ref_exists,
        unsigned length_descriptor, unsigned num_elements, RemappedVariant& remapper_variant) = 0;
    virtual bool get_valid_median(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, void* output_ptr) = 0; 
    virtual bool get_valid_sum(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, void* output_ptr) = 0; 
    virtual bool compute_valid_element_wise_sum(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, void* output_ptr, unsigned num_elements) = 0;
    virtual bool collect_and_extend_fields(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, const void ** output_ptr, unsigned& num_elements) = 0;
};

//Big bag handler functions useful for handling different types of fields (int, char etc)
//Helps avoid writing a whole bunch of switch statements
template<class DataType>
class VariantFieldHandler : public VariantFieldHandlerBase
{
  public:
    VariantFieldHandler() : VariantFieldHandlerBase()
    { 
      //resize once, re-use many times - avoid reallocs()
      m_num_calls_with_valid_data.resize(100u);
      m_bcf_missing_value = get_bcf_missing_value<DataType>();
      //Vector to hold data values for computing median - avoid frequent re-allocs
      m_median_compute_vector.resize(100u);
    }
    ~VariantFieldHandler() = default;
    /*
     * Wrapper function to remap order of elements in fields which depend on order of alleles
     * E.g. PL, AD etc
     */
    virtual void remap_vector_data(std::unique_ptr<VariantFieldBase>& orig_field_ptr, uint64_t curr_call_idx_in_variant, 
        const CombineAllelesLUT& alleles_LUT, unsigned num_merged_alleles, bool non_ref_exists,
        unsigned length_descriptor, unsigned num_merged_elements, RemappedVariant& remapper_variant);
    /*
     * Computes median for a given field over all Calls (only considers calls with valid field)
     */
    virtual bool get_valid_median(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, void* output_ptr); 
    /*
     * Computes sum for a given field over all Calls (only considers calls with valid field)
     */
    virtual bool get_valid_sum(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, void* output_ptr);
    /*
     * Computes element-wise sum for a given field over all Calls (only considers calls with valid field)
     */
    virtual bool compute_valid_element_wise_sum(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, void* output_ptr, unsigned num_elements);
    /*
     * Create an extended vector for use in BCF format fields, return result in output_ptr and num_elements
     */
    bool collect_and_extend_fields(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, const void ** output_ptr, unsigned& num_elements);
  private:
    std::vector<uint64_t> m_num_calls_with_valid_data;
    DataType m_bcf_missing_value;
    //Vector to hold data values for computing median - avoid frequent re-allocs
    std::vector<DataType> m_median_compute_vector;
    //Vector to hold extended vector to use in BCF format fields
    std::vector<DataType> m_extended_field_vector;
};

/*
 * Copies info in Variant object into its result vector
 */
class GA4GHOperator : public SingleVariantOperatorBase
{
  public:
    GA4GHOperator(const VariantQueryConfig& query_config);
    virtual void operate(Variant& variant, const VariantQueryConfig& query_config);
    const Variant& get_remapped_variant() const { return m_remapped_variant; }
    Variant& get_remapped_variant() { return m_remapped_variant; }
    void copy_back_remapped_fields(Variant& variant) const;
  protected:
    Variant m_remapped_variant;
    //Query idxs of fields that need to be remmaped - PL, AD etc
    std::vector<unsigned> m_remapped_fields_query_idxs;
    //Query idx of GT field, could be UNDEFINED_ATTRIBUTE_IDX_VALUE
    unsigned m_GT_query_idx;
    //Get handler based on type of field
    std::unique_ptr<VariantFieldHandlerBase>& get_handler_for_type(std::type_index ty);
    //Handlers for various fields
    std::vector<std::unique_ptr<VariantFieldHandlerBase>> m_field_handlers;
};

class SingleCellOperatorBase
{
  public:
    SingleCellOperatorBase() { ; }
    virtual void operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema)  { ; }
};

class ColumnHistogramOperator : public SingleCellOperatorBase
{
  public:
    ColumnHistogramOperator(uint64_t begin, uint64_t end, uint64_t bin_size);
    virtual void operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema);
    bool equi_partition_and_print_bins(uint64_t num_bins, std::ostream& fptr=std::cout) const; 
  private:
    std::vector<uint64_t> m_bin_counts_vector;
    uint64_t m_begin_column;
    uint64_t m_end_column;
    uint64_t m_bin_size;
};

class VariantCallPrintOperator : public SingleCellOperatorBase
{
  public:
    VariantCallPrintOperator(std::ostream& fptr=std::cout, const std::string& indent_prefix="")
      : SingleCellOperatorBase(), m_fptr(&fptr), m_indent_prefix(indent_prefix)
      {
        m_num_calls_printed = 0ull;
      }
    virtual void operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema)
    {
      if(m_num_calls_printed > 0ull)
        (*m_fptr) << ",\n";
      call.print(*m_fptr, &query_config, m_indent_prefix);
      ++m_num_calls_printed;
    }
  private:
    uint64_t m_num_calls_printed;
    std::string m_indent_prefix;
    std::ostream* m_fptr;
};

//Dump CSV
class VariantCallPrintCSVOperator : public SingleCellOperatorBase
{
  public:
    VariantCallPrintCSVOperator(std::ostream& fptr=std::cout)
      : SingleCellOperatorBase(), m_fptr(&fptr)
      {
      }
    virtual void operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema);
  private:
    std::ostream* m_fptr;
};

/*
 * If the call's column is before the current_start_position, then REF is not valid, set it to "N" (unknown/don't care)
 */
void modify_reference_if_in_middle(VariantCall& curr_call, const VariantQueryConfig& query_config, uint64_t current_start_position);

#endif
