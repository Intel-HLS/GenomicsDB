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

#include <gtest/gtest.h>
#include "variant_operations.h"

class RemapDataTestClass : public RemappedDataWrapperBase
{
   public:
     RemapDataTestClass(const std::vector<std::string>& golden_genotype_combinations)
       : m_golden_genotype_combinations(golden_genotype_combinations),
       m_golden_input_genotype_combinations(golden_genotype_combinations)
     {
     }
     void* put_address(uint64_t input_call_idx, unsigned allele_or_gt_idx)
     {
       throw VariantOperationException("Unimplemented");
     }
     const std::string& get_genotype_combination(const unsigned idx) const
     {
       assert(idx < m_golden_genotype_combinations.size());
       return m_golden_genotype_combinations[idx];
     }
     const std::string& get_input_genotype_combination(const unsigned idx) const
     {
       assert(idx < m_golden_input_genotype_combinations.size());
       return m_golden_input_genotype_combinations[idx];
     }
     void set_golden_input_genotype_combinations(const std::vector<std::string>&  vec)
     {
       m_golden_input_genotype_combinations = vec;
     }
   private:
     std::vector<std::string> m_golden_genotype_combinations;
     std::vector<std::string> m_golden_input_genotype_combinations;
};


template<class DataType>
void verify_remapped_gt_order(const std::vector<DataType>& input_data,
    const uint64_t input_call_idx,
    const CombineAllelesLUT& alleles_LUT,
    const unsigned num_merged_alleles, bool NON_REF_exists, 
    const bool curr_genotype_combination_contains_missing_allele_for_input,
    const unsigned ploidy,
    RemappedDataWrapperBase& remapped_data,
    std::vector<uint64_t>& num_calls_with_valid_data, DataType missing_value,
    const std::vector<int>& remapped_allele_idx_vec_for_current_gt_combination,
    const uint64_t remapped_gt_idx,
    std::vector<int>& input_call_allele_idx_vec
    )
{
  std::string s;
  for(auto i=0u;i<ploidy;++i)
    s += (static_cast<char>(remapped_allele_idx_vec_for_current_gt_combination[i]) + 'A');
  EXPECT_EQ(dynamic_cast<RemapDataTestClass&>(remapped_data).get_genotype_combination(remapped_gt_idx),
      s);
  s.clear();
  auto input_gt_idx = VariantOperations::get_genotype_index(input_call_allele_idx_vec, false);
  auto found_unmatched_allele = false;
  for(auto i=0u;i<ploidy;++i)
  {
    if(input_call_allele_idx_vec[i] < 0)
    {
      EXPECT_TRUE(curr_genotype_combination_contains_missing_allele_for_input);
      found_unmatched_allele = true;
    }
    else
      s += (static_cast<char>(input_call_allele_idx_vec[i]) + 'A');
  }
  if(!curr_genotype_combination_contains_missing_allele_for_input)
  {
    //Verifies functionality of get_genotype_index()
    EXPECT_EQ(dynamic_cast<RemapDataTestClass&>(remapped_data).get_genotype_combination(input_gt_idx),
        s);
    //Verifies that lookup from merged genotype to input genotype is correct
    EXPECT_EQ(dynamic_cast<RemapDataTestClass&>(remapped_data).get_input_genotype_combination(remapped_gt_idx),
        s);
  }
  else
  {
    EXPECT_TRUE(found_unmatched_allele);
    EXPECT_EQ(dynamic_cast<RemapDataTestClass&>(remapped_data).get_input_genotype_combination(remapped_gt_idx),
        ".");
  }
}

void initialize_LUT_identical(CombineAllelesLUT& lut)
{
  lut.resize_luts_if_needed(1, 5);
  lut.reset_luts();
  lut.add_input_merged_idx_pair(0, 0, 0);
  lut.add_input_merged_idx_pair(0, 1, 1);
  lut.add_input_merged_idx_pair(0, 2, 2);
  lut.add_input_merged_idx_pair(0, 3, 3);
  lut.add_input_merged_idx_pair(0, 4, 4);
}

void initialize_LUT_reordered(CombineAllelesLUT& lut)
{
  //Input call has 3 alleles - A,B,<NON_REF>
  lut.resize_luts_if_needed(1, 5);
  lut.reset_luts();
  //A_inp mapped to A_merged
  lut.add_input_merged_idx_pair(0, 0, 0);
  //B_inp mapped to C_merged
  lut.add_input_merged_idx_pair(0, 1, 2);
  //C_inp [NON_REF] mapped to D_merged
  lut.add_input_merged_idx_pair(0, 2, 3); //NON_REF
  //Implies that any alleles corresponding to B_merged will get mapped to NON_REF allele
}

void initialize_LUT_reordered_without_NON_REF(CombineAllelesLUT& lut)
{
  //Input call has 3 alleles - A,B,<NON_REF>
  lut.resize_luts_if_needed(1, 5);
  lut.reset_luts();
  //A_inp mapped to A_merged
  lut.add_input_merged_idx_pair(0, 0, 0);
  //B_inp mapped to C_merged
  lut.add_input_merged_idx_pair(0, 1, 2);
}

//Golden outputs copied from:
//http://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes

TEST(genotype_ordering, haploid_5_alleles)
{
  std::vector<int> tmp_vector;
  std::vector<std::pair<int, int> > tmp_stack;
  std::vector<int> tmp_input_call_vector;
  CombineAllelesLUT tmp_lut;
  initialize_LUT_identical(tmp_lut);
  RemapDataTestClass test_data_obj({"A", "B", "C", "D", "E"});
  std::vector<uint64_t> tmp_count_vector;
  VariantOperations::remap_data_based_on_genotype_general<int>(tmp_vector,
      0ull,
      tmp_lut,
      5, false, 1u,
      test_data_obj,
      tmp_count_vector, bcf_int32_missing,
      tmp_vector, tmp_stack,
      tmp_input_call_vector,
      verify_remapped_gt_order);
}

TEST(genotype_ordering, diploid_2_alleles)
{
  std::vector<int> tmp_vector;
  std::vector<std::pair<int, int> > tmp_stack;
  std::vector<int> tmp_input_call_vector;
  CombineAllelesLUT tmp_lut;
  initialize_LUT_identical(tmp_lut);
  RemapDataTestClass test_data_obj({
      "AA",
      "AB",
      "BB"
      });
  std::vector<uint64_t> tmp_count_vector;
  VariantOperations::remap_data_based_on_genotype_general<int>(tmp_vector,
      0ull,
      tmp_lut,
      2, false, 2u,
      test_data_obj,
      tmp_count_vector, bcf_int32_missing,
      tmp_vector, tmp_stack,
      tmp_input_call_vector,
      verify_remapped_gt_order);
}

TEST(genotype_ordering, diploid_3_alleles)
{
  std::vector<int> tmp_vector;
  std::vector<std::pair<int, int> > tmp_stack;
  std::vector<int> tmp_input_call_vector;
  CombineAllelesLUT tmp_lut;
  initialize_LUT_identical(tmp_lut);
  RemapDataTestClass test_data_obj({
      "AA",
      "AB",
      "BB",
      "AC",
      "BC",
      "CC"
      });
  std::vector<uint64_t> tmp_count_vector;
  VariantOperations::remap_data_based_on_genotype_general<int>(tmp_vector,
      0ull,
      tmp_lut,
      3, false, 2u,
      test_data_obj,
      tmp_count_vector, bcf_int32_missing,
      tmp_vector, tmp_stack,
      tmp_input_call_vector,
      verify_remapped_gt_order);
}

TEST(genotype_ordering, diploid_4_alleles)
{
  std::vector<int> tmp_vector;
  std::vector<std::pair<int, int> > tmp_stack;
  std::vector<int> tmp_input_call_vector;
  CombineAllelesLUT tmp_lut;
  initialize_LUT_identical(tmp_lut);
  RemapDataTestClass test_data_obj({
      "AA",
      "AB",
      "BB",
      "AC",
      "BC",
      "CC",
      "AD",
      "BD",
      "CD",
      "DD"
      });
  std::vector<uint64_t> tmp_count_vector;
  VariantOperations::remap_data_based_on_genotype_general<int>(tmp_vector,
      0ull,
      tmp_lut,
      4, false, 2u,
      test_data_obj,
      tmp_count_vector, bcf_int32_missing,
      tmp_vector, tmp_stack,
      tmp_input_call_vector,
      verify_remapped_gt_order);
}

TEST(genotype_ordering, triploid_2_alleles)
{
  std::vector<int> tmp_vector;
  std::vector<std::pair<int, int> > tmp_stack;
  std::vector<int> tmp_input_call_vector;
  CombineAllelesLUT tmp_lut;
  initialize_LUT_identical(tmp_lut);
  RemapDataTestClass test_data_obj({
      "AAA",
      "AAB",
      "ABB",
      "BBB"
      });
  std::vector<uint64_t> tmp_count_vector;
  VariantOperations::remap_data_based_on_genotype_general<int>(tmp_vector,
      0ull,
      tmp_lut,
      2, false, 3u,
      test_data_obj,
      tmp_count_vector, bcf_int32_missing,
      tmp_vector, tmp_stack,
      tmp_input_call_vector,
      verify_remapped_gt_order);
}

TEST(genotype_ordering, triploid_3_alleles)
{
  std::vector<int> tmp_vector;
  std::vector<std::pair<int, int> > tmp_stack;
  std::vector<int> tmp_input_call_vector;
  CombineAllelesLUT tmp_lut;
  initialize_LUT_identical(tmp_lut);
  RemapDataTestClass test_data_obj({
      "AAA",
      "AAB",
      "ABB",
      "BBB",
      "AAC",
      "ABC",
      "BBC",
      "ACC",
      "BCC",
      "CCC"
      });
  std::vector<uint64_t> tmp_count_vector;
  VariantOperations::remap_data_based_on_genotype_general<int>(tmp_vector,
      0ull,
      tmp_lut,
      3, false, 3u,
      test_data_obj,
      tmp_count_vector, bcf_int32_missing,
      tmp_vector, tmp_stack,
      tmp_input_call_vector,
      verify_remapped_gt_order);
}

TEST(genotype_ordering, triploid_4_alleles)
{
  std::vector<int> tmp_vector;
  std::vector<std::pair<int, int> > tmp_stack;
  std::vector<int> tmp_input_call_vector;
  CombineAllelesLUT tmp_lut;
  initialize_LUT_identical(tmp_lut);
  RemapDataTestClass test_data_obj({
      "AAA",
      "AAB",
      "ABB",
      "BBB",
      "AAC",
      "ABC",
      "BBC",
      "ACC",
      "BCC",
      "CCC",
      "AAD",
      "ABD",
      "BBD",
      "ACD",
      "BCD",
      "CCD",
      "ADD",
      "BDD",
      "CDD",
      "DDD"
      });
  std::vector<uint64_t> tmp_count_vector;
  VariantOperations::remap_data_based_on_genotype_general<int>(tmp_vector,
      0ull,
      tmp_lut,
      4, false, 3u,
      test_data_obj,
      tmp_count_vector, bcf_int32_missing,
      tmp_vector, tmp_stack,
      tmp_input_call_vector,
      verify_remapped_gt_order);
}

TEST(genotype_ordering, triploid_4_alleles_reordered_alleles)
{
  std::vector<int> tmp_vector;
  std::vector<std::pair<int, int> > tmp_stack;
  std::vector<int> tmp_input_call_vector;
  CombineAllelesLUT tmp_lut;
  initialize_LUT_reordered(tmp_lut);
  RemapDataTestClass test_data_obj({
      "AAA",
      "AAB",
      "ABB",
      "BBB",
      "AAC",
      "ABC",
      "BBC",
      "ACC",
      "BCC",
      "CCC",
      "AAD",
      "ABD",
      "BBD",
      "ACD",
      "BCD",
      "CCD",
      "ADD",
      "BDD",
      "CDD",
      "DDD"
      });
  ////A_inp mapped to A_merged
  ////B_inp mapped to C_merged
  ////C_inp [NON_REF] mapped to D_merged
  ////Implies that any alleles corresponding to B_merged will get mapped to NON_REF allele
  test_data_obj.set_golden_input_genotype_combinations({
      "AAA",
      "AAC", //AAB
      "ACC", //ABB
      "CCC", //BBB
      "AAB", //AAC
      "ABC", //ABC -> ACB -> ABC
      "BCC", //BBC
      "ABB", //ACC
      "BBC", //BCC
      "BBB", //CCC
      "AAC", //AAD
      "ACC", //ABD
      "CCC", //BBD
      "ABC", //ACD
      "BCC", //BCD -> CBC -> BCC
      "BBC", //CCD ->
      "ACC", //ADD
      "CCC", //BDD
      "BCC", //CDD
      "CCC"  //DDD
      });
  std::vector<uint64_t> tmp_count_vector;
  VariantOperations::remap_data_based_on_genotype_general<int>(tmp_vector,
      0ull,
      tmp_lut,
      4, true, 3u,
      test_data_obj,
      tmp_count_vector, bcf_int32_missing,
      tmp_vector, tmp_stack,
      tmp_input_call_vector,
      verify_remapped_gt_order);
}

TEST(genotype_ordering, triploid_4_alleles_reordered_alleles_without_NON_REF)
{
  std::vector<int> tmp_vector;
  std::vector<std::pair<int, int> > tmp_stack;
  std::vector<int> tmp_input_call_vector;
  CombineAllelesLUT tmp_lut;
  initialize_LUT_reordered_without_NON_REF(tmp_lut);
  RemapDataTestClass test_data_obj({
      "AAA",
      "AAB",
      "ABB",
      "BBB",
      "AAC",
      "ABC",
      "BBC",
      "ACC",
      "BCC",
      "CCC",
      "AAD",
      "ABD",
      "BBD",
      "ACD",
      "BCD",
      "CCD",
      "ADD",
      "BDD",
      "CDD",
      "DDD"
      });
  //Input call has 2 alleles - A,B
  //A_inp mapped to A_merged
  //B_inp mapped to C_merged
  test_data_obj.set_golden_input_genotype_combinations({
      "AAA",
      ".",   //AAB
      ".",   //ABB
      ".",   //BBB
      "AAB", //AAC
      ".",   //ABC
      ".",   //BBC
      "ABB", //ACC
      ".",   //BCC
      "BBB", //CCC
      ".",   //AAD
      ".", //ABD
      ".", //BBD
      ".", //ACD
      ".", //BCD
      ".", //CCD
      ".", //ADD
      ".", //BDD
      ".", //CDD
      "."  //DDD
      });
  std::vector<uint64_t> tmp_count_vector;
  VariantOperations::remap_data_based_on_genotype_general<int>(tmp_vector,
      0ull,
      tmp_lut,
      4, true, 3u,
      test_data_obj,
      tmp_count_vector, bcf_int32_missing,
      tmp_vector, tmp_stack,
      tmp_input_call_vector,
      verify_remapped_gt_order);
}
