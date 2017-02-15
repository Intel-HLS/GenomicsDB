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

#ifndef GT_LUT_H
#define GT_LUT_H

#include <assert.h>
#include <vector>
#include <stdlib.h>

#define lut_missing_value -1ll
/**
 * LUT = Look Up Table (to avoid confusion with map, unordered_map etc)
 * @brief Base class to store look up information between fields of merged header and input headers
 * @note This is the helper class for VariantHeaderMerger to store mapping for fields and samples
 * Each LUTBase object contains 2 matrices (vector of vector): one for mapping input field idx to merged field idx (m_inputs_2_merged_lut)
 * and the second for mapping merged field idx to input field idx (m_merged_2_inputs_lut).
 *
 * Missing field information is stored as lut_missing_value, but should be checked with is_missing_value() function
 * 
 * The boolean template parameters specify how the 2 tables are laid out in memory - whether the outer vector corresponds to fields or input vcfs.
 * For example, in object of type LUTBase<true, true>, both LUTs are laid out such that m_inputs_2_merged_lut[0] contains mappings 
 * for all fields for input VCF file 0. This would lead to fast traversal of all fields for a given input VCF (cache locality). 
 * However, traversing over all input VCFs for a given field would be slow (many cache misses).
 * The object LUTBase<false,false> would have the exact opposite behavior
 * 
 * The 'best' value of the template parameters depends on the application using the LUT.
 * Almost all the 'complexity' of the code comes from being able to handle the different layouts in a transparent manner
 *
 * Alternate explanation:
 * This class contains two matrices (vector<vector<int64_t>>) to store the mapping information:
 * m_inputs_2_merged_lut and m_merged_2_inputs_lut. You can layout each matrix in one of the 2 following ways:
 * (a) LUT[i][j]  corresponds to input VCF i and field j 
 * (b) LUT[i][j]  corresponds to field i and input VCF j
 * Option (a) is optimal where you are looking at all the fields of a VCF in quick succession,
 * while (b) is optimal when you are looking at all VCFs for a particular field.
 * The 2 boolean template parameters control the layout of the two matrices. If the parameter value is true,
 * then option (a) is picked, else option (b)
 *
 * Although the class provides functions to resize the tables, for obtaining good performance, reallocations should be extremely
 * infrequent. Making the resize_luts_if_needed() a protected member forces developers to think twice instead of blindly calling this function.
 *
 * Uses the enable_if trick from http://en.cppreference.com/w/cpp/types/enable_if  (foo3) to handle different memory layouts
 * 
 **/
template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
class LUTBase
{
  public:
  /**
   * @brief: clear all mappings
   */
  inline void reset_luts()
  {
    for(auto& vec : m_inputs_2_merged_lut)
      reset_vector(vec);
    for(auto& vec : m_merged_2_inputs_lut)
      reset_vector(vec);
  }

  /*
   * @brief Add a valid mapping between input VCF and merged VCF
   * @note all parameters should be valid parameters, no lut_missing_value, use reset_() functions to invalidate existing mapping
   * @param inputGVCFIdx index of the input VCF file
   * @param inputIdx index of the field in the input VCF file - field could be anything header field,sample,allele etc
   * @param mergedIdx index of the field in the merged VCF file
   */
  inline void add_input_merged_idx_pair(int64_t inputGVCFIdx, int64_t inputIdx, int64_t mergedIdx)
  {
    set_merged_idx_for_input(inputGVCFIdx, inputIdx, mergedIdx);
    set_input_idx_for_merged(inputGVCFIdx, inputIdx, mergedIdx);
  }

  /**
   * @brief Get field idx for input VCF inputGVCFIdx corresponding to field idx mergedIdx in the mergedVCF file
   * @note Uses the enable_if trick from http://en.cppreference.com/w/cpp/types/enable_if  (foo3) to handle different memory layouts
   * The enable_if<M> corresponds to the case where merged_2_inputs_LUT_is_input_ordered = true, hence, the rows correspond to input VCFs
   * The enable_if<!M> corresponds to the case where merged_2_inputs_LUT_is_input_ordered = false, hence, the rows correspond to fields
   * @param inputGVCFIdx index of the input VCF file
   * @param mergedIdx index of the field in the merged VCF file
   * @return index of the field in the input VCF file
   */
  template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
  inline int64_t get_input_idx_for_merged(int64_t inputGVCFIdx, int64_t mergedIdx) const
  { return get_lut_value(m_merged_2_inputs_lut, inputGVCFIdx, mergedIdx); }
  template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
  inline int64_t get_input_idx_for_merged(int64_t inputGVCFIdx, int64_t mergedIdx) const
  { return get_lut_value(m_merged_2_inputs_lut, mergedIdx, inputGVCFIdx); }

  /**
   * @brief Get field idx for the merged VCF corresponding to field idx inputIdx in the input VCF of index inputGVCFIdx
   * @note Uses the enable_if trick from http://en.cppreference.com/w/cpp/types/enable_if  (foo3) to handle different memory layouts
   * The enable_if<M> corresponds to the case where inputs_2_merged_LUT_is_input_ordered = true, hence, the rows correspond to input VCFs
   * The enable_if<!M> corresponds to the case where inputs_2_merged_LUT_is_input_ordered = false, hence, the rows correspond to fields
   * @param inputGVCFIdx index of the input VCF file
   * @param inputIdx index of the field in the input VCF file
   * @return index of the field in the merged VCF file
   */
  template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
  inline int64_t get_merged_idx_for_input(int64_t inputGVCFIdx, int64_t inputIdx) const
  { return get_lut_value(m_inputs_2_merged_lut, inputGVCFIdx, inputIdx); }
  template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
  inline int64_t get_merged_idx_for_input(int64_t inputGVCFIdx, int64_t inputIdx) const
  { return get_lut_value(m_inputs_2_merged_lut, inputIdx, inputGVCFIdx); }

  /**
   * @brief reset/invalidate merged field index for field inputIdx of input VCF inputGVCFIdx
   * @param inputGVCFIdx index of the input VCF file
   * @param inputIdx index of the field in the input VCF file
   */
  inline void reset_merged_idx_for_input(int64_t inputGVCFIdx, int64_t inputIdx)
  { set_merged_idx_for_input(inputGVCFIdx, inputIdx, lut_missing_value); }
  /**
   * @brief reset/invalidate the input field index for input VCF inputGVCFIdx for merged field mergedIdx
   * @param inputGVCFIdx index of the input VCF file
   * @param mergedIdx index of the field in the merged VCF file
   */
  inline void reset_input_idx_for_merged(int64_t inputGVCFIdx, int64_t mergedIdx)
  { set_input_idx_for_merged(inputGVCFIdx, lut_missing_value, mergedIdx); }

  static inline bool is_missing_value(int64_t value) { return value == lut_missing_value; }

  protected:
  //Only inherited classes should call constructor,destructor etc
  LUTBase(); 
  LUTBase(int64_t numInputGVCFs, int64_t numMergedFields);
  ~LUTBase() = default;
  /**
   * @brief deallocates memory
   */
  void clear();

  int64_t m_num_input_vcfs;
  int64_t m_num_merged_fields;

  /**
   *  @brief resize LUT functions 
   *  @note should be called relatively infrequently (more precisely, the reallocation code inside these resize functions should be called
   *  infrequently
   *  @param numInputGVCFs number of input VCFs
   *  @param numMergedFields number of fields combined across all input VCFs
   */
  template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
  void resize_inputs_2_merged_lut_if_needed(int64_t numInputGVCFs, int64_t numMergedFields)
  { resize_and_reset_lut(m_inputs_2_merged_lut, numInputGVCFs, numMergedFields, m_num_input_vcfs, m_num_merged_fields); }

  template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
  void resize_inputs_2_merged_lut_if_needed(int64_t numInputGVCFs, int64_t numMergedFields)
  { resize_and_reset_lut(m_inputs_2_merged_lut, numMergedFields, numInputGVCFs, m_num_merged_fields, m_num_input_vcfs); }

  template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
  void resize_merged_2_inputs_lut_if_needed(int64_t numInputGVCFs, int64_t numMergedFields)
  { resize_and_reset_lut(m_merged_2_inputs_lut, numInputGVCFs, numMergedFields, m_num_input_vcfs, m_num_merged_fields); }

  template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
  void resize_merged_2_inputs_lut_if_needed(int64_t numInputGVCFs, int64_t numMergedFields)
  { resize_and_reset_lut(m_merged_2_inputs_lut, numMergedFields, numInputGVCFs, m_num_merged_fields, m_num_input_vcfs); }

  /*
   * @brief wrapper around single LUT resize functions
   */
  void resize_luts_if_needed(int64_t numInputGVCFs, int64_t numMergedFields)
  {
    resize_merged_2_inputs_lut_if_needed(numInputGVCFs, numMergedFields);
    resize_inputs_2_merged_lut_if_needed(numInputGVCFs, numMergedFields);
  }
  private:
  //why not unordered_map? because I feel the need, the need for speed
  std::vector<std::vector<int64_t>> m_inputs_2_merged_lut;
  std::vector<std::vector<int64_t>> m_merged_2_inputs_lut;
  /**
   * @brief invalidate/reset all mappings in a vector
   * @note sets all elements to missing
   * @param vec the vector to reset
   * @param from offset in the vector from which to start reset, 0 by default
   */
  void reset_vector(std::vector<int64_t>& vec, int64_t from=0u);
  /**
   * @brief resize and reset a vector
   * @note resize and reset is done only if new_size > vec.size()
   */
  void resize_and_reset_vector(std::vector<int64_t>& vec, int64_t new_size);
  /**
   * @brief resize and reset a LUT
   * @note resize and reset is done only if new_size > old_size
   */
  void resize_and_reset_lut(std::vector<std::vector<int64_t>>& lut, int64_t new_lut_size, int64_t new_size, int64_t& numRowsVar, int64_t& numColsVar);

  /**
   * @brief get LUT value at a particular row,column
   * @note should be called only from the public wrapper functions get_*() as the wrappers take care of memory layout
   * @param lut LUT to access
   * @param rowIdx row
   * @param columnIdx column
   * @return value at lut[row][column], could be invalid, check with is_missing()
   */
  inline int64_t get_lut_value(const std::vector<std::vector<int64_t>>& lut, int64_t rowIdx, int64_t columnIdx) const
  {
    assert(rowIdx >= 0);
    assert(rowIdx < static_cast<int64_t>(lut.size()));
    assert(columnIdx >= 0);
    assert(columnIdx < static_cast<int64_t>(lut[rowIdx].size()));
    return lut[rowIdx][columnIdx];
  }

  /**
   * @brief set LUT value at a particular row,column
   * @note should be called only from the public wrapper functions add_input_merged_idx_pair() or reset_*() as the wrappers take care of memory layout
   * @param lut LUT to access
   * @param rowIdx row
   * @param columnIdx column
   * @param value value to write at lut[row][column] 
   */
  inline void set_lut_value(std::vector<std::vector<int64_t>>& lut, int64_t rowIdx, int64_t columnIdx, int64_t value)
  {
    assert(rowIdx >= 0);
    assert(rowIdx < static_cast<int64_t>(lut.size()));
    assert(columnIdx >= 0);
    assert(columnIdx < static_cast<int64_t>(lut[rowIdx].size()));
    lut[rowIdx][columnIdx] = value;
  }

  /**
   * @brief set merged field idx value (mergedIdx) corresponding to field idx inputIdx for input VCF inputGVCFIdx
   * @note should be called only from the public wrapper functions add_input_merged_idx_pair() or reset_*() as the wrappers take care of memory layout
   * @param inputGVCFIdx index of the input VCF file
   * @param inputIdx index of the field in the input VCF file
   * @param mergedIdx index of the field in the merged VCF file
   */
  template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
  inline void set_merged_idx_for_input(int64_t inputGVCFIdx, int64_t inputIdx, int64_t mergedIdx)
  { set_lut_value(m_inputs_2_merged_lut, inputGVCFIdx, inputIdx, mergedIdx); } 

  template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
  inline void set_merged_idx_for_input(int64_t inputGVCFIdx, int64_t inputIdx, int64_t mergedIdx)
  { set_lut_value(m_inputs_2_merged_lut, inputIdx, inputGVCFIdx, mergedIdx); } 

  /**
   * @brief set input field idx value (inputIdx) for input VCF inputGVCFIdx corresponding to field idx mergedIdx in the merged VCF
   * @note should be called only from the public wrapper functions add_input_merged_idx_pair() or reset_*() as the wrappers take care of memory layout
   * @param inputGVCFIdx index of the input VCF file
   * @param inputIdx index of the field in the input VCF file
   * @param mergedIdx index of the field in the merged VCF file
   */
  template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
  inline void set_input_idx_for_merged(int64_t inputGVCFIdx, int64_t inputIdx, int64_t mergedIdx)
  { set_lut_value(m_merged_2_inputs_lut, inputGVCFIdx, mergedIdx, inputIdx); }

  template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
  inline void set_input_idx_for_merged(int64_t inputGVCFIdx, int64_t inputIdx, int64_t mergedIdx)
  { set_lut_value(m_merged_2_inputs_lut, mergedIdx, inputGVCFIdx, inputIdx); }

};


template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
class GoldLUTTemplate
: public LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>
{
  public:
    GoldLUTTemplate(int64_t numInputGVCFs, int64_t numEntries)
      : LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>(numInputGVCFs, numEntries)
    { }
    void clear()
    { LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::clear(); }
    inline void resize_luts_if_needed(int64_t numInputGVCFs, int64_t numEntries)
    {
      LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::resize_luts_if_needed(
          numInputGVCFs, numEntries); 
    }
    inline int64_t get_gold_idx_for_test(int64_t inputGVCFIdx, int64_t inputIdx)
    { return LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::get_merged_idx_for_input(inputGVCFIdx, inputIdx); }
    inline int64_t get_test_idx_for_gold(int64_t inputGVCFIdx, int64_t goldIdx)
    { return LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::get_input_idx_for_merged(inputGVCFIdx, goldIdx); }
};

using GoldLUT = GoldLUTTemplate<true, true>;

/**
 * @brief LUT class for storing mappings between allele vectors in the merged file and input VCF files
 * Since the #alleles per site is expected to be small, this class sets the number of fields to 10. This makes any subsequent re-allocations
 * unlikely. The function resize_luts_if_needed() will almost always return immediately after failing the if condition
 */
template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
class MergedAllelesIdxLUT
: public LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>
{
  private:
    static const auto m_DEFAULT_NUM_INPUT_GVCFS=10u;
    static const auto m_DEFAULT_INIT_NUM_ALLELES=10u;
  public:
    MergedAllelesIdxLUT()
      : LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>(m_DEFAULT_NUM_INPUT_GVCFS,
          m_DEFAULT_INIT_NUM_ALLELES)
    { m_max_num_alleles = m_DEFAULT_INIT_NUM_ALLELES; }
    MergedAllelesIdxLUT(int64_t numInputGVCFs)
      : LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>(numInputGVCFs,
          m_DEFAULT_INIT_NUM_ALLELES)
    { m_max_num_alleles = m_DEFAULT_INIT_NUM_ALLELES; }
    inline void resize_luts_if_needed(int64_t numMergedAlleles)
    {
      if(numMergedAlleles > m_max_num_alleles)
      {
        LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::resize_luts_if_needed(
            LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::m_num_input_vcfs, numMergedAlleles); 
        m_max_num_alleles = numMergedAlleles;
      }
    }
    inline void resize_luts_if_needed(int64_t numInputGVCFs, int64_t numMergedAlleles)
    {
      LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::resize_luts_if_needed(
          numInputGVCFs, numMergedAlleles); 
    }
    inline static bool is_missing_value(unsigned val)
    {
      return LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::is_missing_value(
          static_cast<int>(val));
    }
  private:
    int64_t m_max_num_alleles;
};

/*NOTE: Needs explicit instantiation in .cpp file to use this type alias*/
using CombineAllelesLUT = MergedAllelesIdxLUT<true,true>;

//Not a full LUT, each matrix is actually just a vector (single row matrix)
class QueryIdxToKnownVariantFieldsEnumLUT : public LUTBase<true, true>
{
  public:
    QueryIdxToKnownVariantFieldsEnumLUT()
      : LUTBase<true, true>(1u, 100u)
    {
      m_LUT_size = 100u;
    }
    inline void resize_luts_if_needed(int64_t numQueriedAttributes, int64_t numKnownVariantFields)
    {
      int64_t maxValue = std::max(numKnownVariantFields, numQueriedAttributes);
      if(maxValue > m_LUT_size)
      {
        LUTBase<true, true>::resize_luts_if_needed(1u, maxValue); 
        m_LUT_size = maxValue;
      }
    }
    void add_query_idx_known_field_enum_mapping(int64_t queryIdx, int64_t knownFieldEnum)
    {
      add_input_merged_idx_pair(0u, queryIdx, knownFieldEnum);
    }
    unsigned get_known_field_enum_for_query_idx(int64_t queryIdx) const
    {
      return static_cast<unsigned>(get_merged_idx_for_input(0u, queryIdx));
    }
    unsigned get_query_idx_for_known_field_enum(int64_t knownFieldEnum) const
    {
      return static_cast<unsigned>(get_input_idx_for_merged(0u, knownFieldEnum));
    }
    bool is_defined_value(unsigned val) const
    {
      return !is_missing_value(static_cast<int>(val));
    }
  private:
    int64_t m_LUT_size;
};

//Not a full LUT, each matrix is actually just a vector (single row matrix)
class SchemaIdxToKnownVariantFieldsEnumLUT : public LUTBase<true, true>
{
  public:
    SchemaIdxToKnownVariantFieldsEnumLUT()
      : LUTBase<true, true>(1u, 100u)
    {
      m_LUT_size = 100u;
    }
    inline void resize_luts_if_needed(int64_t numSchemaAttributes, int64_t numKnownVariantFields)
    {
      int64_t maxValue = std::max(numKnownVariantFields, numSchemaAttributes);
      if(maxValue > m_LUT_size)
      {
        LUTBase<true, true>::resize_luts_if_needed(1u, maxValue); 
        m_LUT_size = maxValue;
      }
    }
    void add_schema_idx_known_field_mapping(int64_t schemaIdx, int64_t knownFieldEnum)
    {
      add_input_merged_idx_pair(0u, schemaIdx, knownFieldEnum);
    }
    unsigned get_known_field_enum_for_schema_idx(int64_t schemaIdx) const
    {
      return static_cast<unsigned>(get_merged_idx_for_input(0u, schemaIdx));
    }
    unsigned get_schema_idx_for_known_field_enum(int64_t knownFieldEnum) const
    {
      return static_cast<unsigned>(get_input_idx_for_merged(0u, knownFieldEnum));
    }
    bool is_defined_value(unsigned val) const
    {
      return !is_missing_value(static_cast<int>(val));
    }
  private:
    int64_t m_LUT_size;
};

//Not a full LUT, each matrix is actually just a vector (single row matrix)
class QueryIdxToVCFHeaderFieldIdxLUT : public LUTBase<true, true>
{
  public:
    QueryIdxToVCFHeaderFieldIdxLUT()
      : LUTBase<true, true>(1u, 100u)
    {
      m_LUT_size = 100u;
    }
    inline void resize_luts_if_needed(int numQueryAttributes, int numHeaderFields)
    {
      int64_t maxValue = std::max(numHeaderFields, numQueryAttributes);
      if(maxValue > m_LUT_size)
      {
        LUTBase<true, true>::resize_luts_if_needed(1u, maxValue); 
        m_LUT_size = maxValue;
      }
    }
    void add_query_idx_header_field_mapping(unsigned queryIdx, int header_field_idx)
    {
      add_input_merged_idx_pair(0u, queryIdx, header_field_idx);
    }
    int get_header_field_idx_for_query_idx(unsigned queryIdx) const
    {
      return static_cast<int>(get_merged_idx_for_input(0u, queryIdx));
    }
    unsigned get_query_idx_for_header_field_idx(int header_field_idx) const
    {
      return static_cast<unsigned>(get_input_idx_for_merged(0u, header_field_idx));
    }
    bool is_defined_value(unsigned val) const
    {
      return !is_missing_value(static_cast<int>(val));
    }
    bool is_defined_value(int val) const
    {
      return !is_missing_value(val);
    }
  private:
    int64_t m_LUT_size;
};
#endif
