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

#ifndef VARIANT_OPERATIONS_H
#define VARIANT_OPERATIONS_H

#include "variant.h"
#include "lut.h"
#include "variant_cell.h"
#include "json_config.h"

class VariantOperationException : public std::exception {
  public:
    VariantOperationException(const std::string m="") : msg_(std::string("Variant operation exception: ")+m) { ; }
    ~VariantOperationException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

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

template<class DataType>
using remap_operator_function_type = void(*)(const std::vector<DataType>& input_data,
    const uint64_t input_call_idx,
    const CombineAllelesLUT& alleles_LUT,
    const unsigned num_merged_alleles, bool NON_REF_exists,
    const bool curr_genotype_combination_contains_missing_allele_for_input,
    const unsigned ploidy,
    RemappedDataWrapperBase& remapped_data,
    std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value,
    const std::vector<int>& remapped_allele_idx_vec_for_current_gt_combination,
    const uint64_t remapped_gt_idx,
    std::vector<int>& input_call_allele_idx_vec_for_current_gt_combination
    );

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
        const CombineAllelesLUT& alleles_LUT, const uint64_t input_call_idx,
        const unsigned num_merged_alleles, const bool has_NON_REF, const FieldLengthDescriptor& length_descriptor);
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
        const CombineAllelesLUT& alleles_LUT,
        const unsigned num_merged_alleles, bool NON_REF_exists, const unsigned ploidy,
        RemappedDataWrapperBase& remapped_data,
        std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value,
        std::vector<int>& remapped_allele_idx_vec_for_current_gt_combination,
        std::vector<std::pair<int, int> >& ploidy_index_allele_index_stack,
        std::vector<int>& input_call_allele_idx_vec_for_current_gt_combination);
    /*
     * Reorders fields whose length and order depend on the number of genotypes (BCF_VL_G)
     * for haploid data
     */
    template<class DataType>
    static void remap_data_based_on_genotype_haploid(const std::vector<DataType>& input_data,
        const uint64_t input_call_idx, 
        const CombineAllelesLUT& alleles_LUT,
        const unsigned num_merged_alleles, bool NON_REF_exists,
        RemappedDataWrapperBase& remapped_data,
        std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value);
    /*
     * Reorders fields whose length and order depend on the number of genotypes (BCF_VL_G)
     * for diploid data
     */
    template<class DataType>
    static void remap_data_based_on_genotype_diploid(const std::vector<DataType>& input_data,
        const uint64_t input_call_idx, 
        const CombineAllelesLUT& alleles_LUT,
        const unsigned num_merged_alleles, bool NON_REF_exists,
        RemappedDataWrapperBase& remapped_data,
        std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value);
    /*
     * Reorders fields whose length and order depend on the number of genotypes (BCF_VL_G)
     * for general ploidy
     */
    template<class DataType>
    static void remap_data_based_on_genotype_general(const std::vector<DataType>& input_data,
        const uint64_t input_call_idx,
        const CombineAllelesLUT& alleles_LUT,
        const unsigned num_merged_alleles, bool NON_REF_exists, const unsigned ploidy,
        RemappedDataWrapperBase& remapped_data,
        std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value,
        std::vector<int>& remapped_allele_idx_vec_for_current_gt_combination,
        std::vector<std::pair<int, int> >& ploidy_index_allele_index_stack,
        std::vector<int>& input_call_allele_idx_vec_for_current_gt_combination,
        remap_operator_function_type<DataType> op
        );

    static void do_dummy_genotyping(Variant& variant, std::ostream& output);

    template<class DataType>
    static void null_remap_genotype_operator(const std::vector<DataType>& input_data,
        const uint64_t input_call_idx,
        const CombineAllelesLUT& alleles_LUT,
        const unsigned num_merged_alleles, bool NON_REF_exists,
        const bool curr_genotype_combination_contains_missing_allele_for_input,
        const unsigned ploidy,
        RemappedDataWrapperBase& remapped_data,
        std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value,
        const std::vector<int>& remapped_allele_idx_vec_for_current_gt_combination,
        const uint64_t remapped_gt_idx,
        std::vector<int>& input_call_allele_idx_vec_for_current_gt_combination
        )
    {
      return;
    }

    template<class DataType>
    static void reorder_field_based_on_genotype_index(const std::vector<DataType>& input_data,
        const uint64_t input_call_idx,
        const CombineAllelesLUT& alleles_LUT,
        const unsigned num_merged_alleles, bool NON_REF_exists,
        const bool curr_genotype_combination_contains_missing_allele_for_input,
        const unsigned ploidy,
        RemappedDataWrapperBase& remapped_data,
        std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value,
        const std::vector<int>& remapped_allele_idx_vec_for_current_gt_combination,
        const uint64_t remapped_gt_idx,
        std::vector<int>& input_call_allele_idx_vec_for_current_gt_combination
        );

    //Combination operation
    static inline uint64_t nCr(const int64_t n, const int64_t r)
    {
      if(n < 0 || r < 0 || n < r)
        return 0u;
      //Not the fastest method - if this is a bottleneck, implement faster algorithm
      auto begin_idx = 0ull;
      auto end_idx = 0ull;
      if(r > n-r)
      {
        begin_idx = r+1u;
        end_idx = n-r;
      }
      else
      {
        begin_idx = n-r+1u;
        end_idx = r;
      }
      auto numerator = 1ull;
      for(auto i=begin_idx;i<=static_cast<uint64_t>(n);++i)
        numerator *= i;
      auto denominator = 1ull;
      for(auto i=1ull;i<=end_idx;++i)
        denominator *= i;
      return numerator/denominator;
    }
    //Return genotype index given a list of allele indexes
    static uint64_t get_genotype_index(std::vector<int>& allele_idx_vec, const bool is_sorted);
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
      m_is_reference_block_only = false;
    }
    virtual ~SingleVariantOperatorBase() { }
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
    //is pure reference block if REF is 1 char, and ALT contains only <NON_REF>
    bool m_is_reference_block_only;
};

/*
 *  Prints interesting positions in the array
 *  A position with column = X is interesting iff
 *  a. a new variant call begins at X OR
 *  b. a variant call ends at X-1
 */
class InterestingLocationsPrinter : public SingleVariantOperatorBase
{
  public:
    InterestingLocationsPrinter(std::ostream& fptr)
      : SingleVariantOperatorBase()
    {
      m_fptr = &fptr;
    }
    virtual void operate(Variant& variant, const VariantQueryConfig& query_config);
  protected:
    //Output stream
    std::ostream* m_fptr;
};

class MaxAllelesCountOperator : public SingleVariantOperatorBase
{
  public:
    class AlleleTracker
    {
      public:
        AlleleTracker(const std::string& ref, const std::vector<std::string>& alt, const int64_t column)
        {
          m_merged_reference_allele = ref;
          m_merged_alt_alleles = alt;
          m_column = column;
        }
        std::string m_merged_reference_allele;
        std::vector<std::string> m_merged_alt_alleles;
        int64_t m_column;
        size_t size() const { return m_merged_alt_alleles.size(); }
    };
    class AlleleTrackerCompare
    {
      public:
        bool operator()(const AlleleTracker& lhs, const AlleleTracker& rhs) const
        {
          return (lhs.size() > rhs.size());
        }
    };
  public:
    MaxAllelesCountOperator(const unsigned top_count=25u)
      : SingleVariantOperatorBase()
    {
      m_top_count = top_count;
      m_curr_count = 0u;
      m_total_lines = 0ull;
    }
    ~MaxAllelesCountOperator()
    {
      std::cerr << "TOTAL "<<m_total_lines<<"\n";
      while(!m_top_alleles_pq.empty())
      {
        auto& top_value = m_top_alleles_pq.top();
        std::cerr << top_value.m_column << "," << top_value.m_merged_reference_allele
          << "," << top_value.size();
        for(const auto& alt_allele : top_value.m_merged_alt_alleles)
          std::cerr << "," << alt_allele;
        std::cerr << "\n";
        m_top_alleles_pq.pop();
      }
    }
    virtual void operate(Variant& variant, const VariantQueryConfig& query_config)
    {
      SingleVariantOperatorBase::operate(variant, query_config);
      ++m_total_lines;
      if(m_curr_count < m_top_count || m_merged_alt_alleles.size() > m_top_alleles_pq.top().size())
      {
        if(m_curr_count < m_top_count)
          ++m_curr_count;
        else
          m_top_alleles_pq.pop();
        m_top_alleles_pq.push(AlleleTracker(m_merged_reference_allele, m_merged_alt_alleles, variant.get_column_begin()));
      }
    }
  protected:
    unsigned m_top_count;
    unsigned m_curr_count;
    std::priority_queue<AlleleTracker, std::vector<AlleleTracker>, AlleleTrackerCompare> m_top_alleles_pq;
    uint64_t m_total_lines;
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
        const CombineAllelesLUT& alleles_LUT,
        unsigned num_merged_alleles, bool non_ref_exists, const unsigned ploidy,
        const FieldLengthDescriptor& length_descriptor, unsigned num_elements, RemappedVariant& remapper_variant) = 0;
    virtual bool get_valid_median(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, void* output_ptr, unsigned& num_valid_elements) = 0;
    virtual bool get_valid_sum(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, void* output_ptr, unsigned& num_valid_elements) = 0;
    virtual bool get_valid_mean(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, void* output_ptr, unsigned& num_valid_elements) = 0;
    virtual bool compute_valid_element_wise_sum(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, const void** output_ptr, unsigned& num_elements) = 0;
    virtual bool concatenate_field(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, const void** output_ptr, unsigned& num_elements) = 0;
    virtual bool collect_and_extend_fields(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, const void ** output_ptr, uint64_t& num_elements,
        const bool use_missing_values_only_not_vector_end=false, const bool use_vector_end_only=false,
        const bool is_GT_field = false) = 0;
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
      //Vector to hold element-wise operations result
      m_element_wise_operations_result.resize(100u);
    }
    ~VariantFieldHandler() = default;
    /*
     * Wrapper function to remap order of elements in fields which depend on order of alleles
     * E.g. PL, AD etc
     */
    virtual void remap_vector_data(std::unique_ptr<VariantFieldBase>& orig_field_ptr, uint64_t curr_call_idx_in_variant, 
        const CombineAllelesLUT& alleles_LUT,
        unsigned num_merged_alleles, bool non_ref_exists, const unsigned ploidy,
        const FieldLengthDescriptor& length_descriptor, unsigned num_merged_elements, RemappedVariant& remapper_variant);
    /*
     * Computes median for a given field over all Calls (only considers calls with valid field)
     */
    virtual bool get_valid_median(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, void* output_ptr, unsigned& num_valid_elements);
    /*
     * Computes sum for a given field over all Calls (only considers calls with valid field)
     */
    virtual bool get_valid_sum(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, void* output_ptr, unsigned& num_valid_elements);
    /*
     * Computes mean for a given field over all Calls (only considers calls with valid field and with valid elements in the field)
     */
    virtual bool get_valid_mean(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, void* output_ptr, unsigned& num_valid_elements);
    /*
     * Concatenates all elements of a field from all Calls in a Variant into a vector
     * Includes both valid and invalid values
     */
    virtual bool concatenate_field(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, const void** output_ptr, unsigned& num_elements);
    /*
     * Computes element-wise sum for a given field over all Calls (only considers calls with valid field)
     */
    virtual bool compute_valid_element_wise_sum(const Variant& variant, const VariantQueryConfig& query_config,
        unsigned query_idx, const void** output_ptr, unsigned& num_elements);
    /*
     * Create an extended vector for use in BCF format fields, return result in output_ptr and num_elements
     */
    bool collect_and_extend_fields(const Variant& variant, const VariantQueryConfig& query_config, 
        unsigned query_idx, const void ** output_ptr, uint64_t& num_elements,
        const bool use_missing_values_only_not_vector_end=false, const bool use_vector_end_only=false,
        const bool is_GT_field=false);
  private:
    std::vector<uint64_t> m_num_calls_with_valid_data;
    DataType m_bcf_missing_value;
    //Vector to hold data values for computing median - avoid frequent re-allocs
    std::vector<DataType> m_median_compute_vector;
    //Vector to hold extended vector to use in BCF format fields
    std::vector<DataType> m_extended_field_vector;
    //Vector to hold data for element wise operations
    std::vector<DataType> m_element_wise_operations_result;
    //Data structures for iterating over genotypes in the non-diploid and non-haploid
    //ploidy. Avoids dynamic memory re-allocation
    //Vector storing allele indexes for current genotype
    std::vector<int> m_allele_idx_vec_for_current_genotype;
    //To avoid dynamic memory allocation - a placeholder for storing
    //input alleles for the given genotype combination in the merged variant
    std::vector<int> m_input_call_allele_idx_vec;
    //"Stack" for enumerating genotypes - rather than using a recursive function
    //as described in http://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes,
    //keep stack of (ploidy index, allele_idx) - only when the stack is empty are all genotypes enumerated
    std::vector<std::pair<int, int> > m_ploidy_index_alleles_index_stack;
};

/*
 * Copies info in Variant object into its result vector
 */
class GA4GHOperator : public SingleVariantOperatorBase
{
  public:
    GA4GHOperator(const VariantQueryConfig& query_config, const unsigned max_diploid_alt_alleles_that_can_be_genotyped=MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED);
    virtual void operate(Variant& variant, const VariantQueryConfig& query_config);
    const Variant& get_remapped_variant() const { return m_remapped_variant; }
    Variant& get_remapped_variant() { return m_remapped_variant; }
    void copy_back_remapped_fields(Variant& variant) const;
    bool too_many_alt_alleles_for_genotype_length_fields(unsigned num_alt_alleles) const { return num_alt_alleles > m_max_diploid_alt_alleles_that_can_be_genotyped; }
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
    //Max alt alleles that can be handled for computing the PL fields - default 50
    unsigned m_max_diploid_alt_alleles_that_can_be_genotyped;
    //1 per VariantCall (row)
    std::vector<unsigned> m_ploidy;
};

class SingleCellOperatorBase
{
  public:
    SingleCellOperatorBase() { ; }
    virtual void operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema)  { ; }
    virtual void operate_on_columnar_cell(const GenomicsDBColumnarCell& cell, const VariantQueryConfig& query_config,
        const VariantArraySchema& schema)
    {
      throw VariantOperationException("Sub-classes should override operate_on_columnar_cell()");
    }
    virtual void finalize() { ; } //do nothing
};

class ColumnHistogramOperator : public SingleCellOperatorBase
{
  public:
    ColumnHistogramOperator(uint64_t begin, uint64_t end, uint64_t bin_size);
    virtual void operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema);
    void operate_on_columnar_cell(const GenomicsDBColumnarCell& cell, const VariantQueryConfig& query_config,
        const VariantArraySchema& schema);
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
    VariantCallPrintOperator(std::ostream& fptr=std::cout, const std::string& indent_prefix="", const VidMapper* vid_mapper=0)
      : SingleCellOperatorBase(), m_fptr(&fptr), m_indent_prefix(indent_prefix)
      {
        m_num_calls_printed = 0ull;
        m_num_query_intervals_printed = 0ull;
        m_vid_mapper = (vid_mapper && vid_mapper->is_initialized()) ? vid_mapper : 0;
        m_indent_prefix_plus_one = m_indent_prefix + g_json_indent_unit;
        m_indent_prefix_plus_two = m_indent_prefix_plus_one + g_json_indent_unit;
      }
    void operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema);
    void operate_on_columnar_cell(const GenomicsDBColumnarCell& cell, const VariantQueryConfig& query_config,
        const VariantArraySchema& schema);
    void finalize();
  private:
    uint64_t m_num_calls_printed;
    uint64_t m_num_query_intervals_printed;
    std::string m_indent_prefix;
    std::string m_indent_prefix_plus_one;
    std::string m_indent_prefix_plus_two;
    std::ostream* m_fptr;
    const VidMapper* m_vid_mapper;
};

//Dump CSV
class VariantCallPrintCSVOperator : public SingleCellOperatorBase
{
  public:
    VariantCallPrintCSVOperator(std::ostream& fptr=std::cout)
      : SingleCellOperatorBase(), m_fptr(&fptr)
      {
      }
    void operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema);
    void operate_on_columnar_cell(const GenomicsDBColumnarCell& cell, const VariantQueryConfig& query_config,
        const VariantArraySchema& schema);
  private:
    std::ostream* m_fptr;
};

class AlleleCountOperator : public SingleCellOperatorBase
{
  public:
    AlleleCountOperator(const VidMapper& vid_mapper, const VariantQueryConfig& query_config);
    void operate(VariantCall& call, const VariantQueryConfig& query_config, const VariantArraySchema& schema);
    void operate_on_columnar_cell(const GenomicsDBColumnarCell& cell, const VariantQueryConfig& query_config,
        const VariantArraySchema& schema);
    void print_allele_counts(std::ostream& fptr=std::cout) const;
    void normalize_REF_ALT_pair(std::pair<std::string, std::string>& REF_ALT_pair);
    std::map<std::pair<std::string, std::string>, uint64_t>& get_REF_ALT_to_count_map(
        const int64_t curr_column);
  private:
    unsigned m_GT_query_idx;
    unsigned m_REF_query_idx;
    unsigned m_ALT_query_idx;
    unsigned m_GT_step_value;
    //Outer vec - 1 per query column position
    //map - 1 per cell begin position
    //map: REF,ALT -> count
    std::vector<std::map<int64_t, std::map<std::pair<std::string, std::string>, uint64_t>>>
      m_column_to_REF_ALT_to_count_vec;
    const VidMapper* m_vid_mapper;
    //Avoid memory allocation
    std::vector<size_t> m_cell_ALT_offsets;
};
/*
 * If the call's column is before the current_start_position, then REF is not valid, set it to "N" (unknown/don't care)
 */
void modify_reference_if_in_middle(VariantCall& curr_call, const VariantQueryConfig& query_config, uint64_t current_start_position);

#endif
