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
#include "vid_mapper.h"
#include "genomicsdb_multid_vector_field.h"

TEST(multid_vector, 2D_test)
{
  FieldInfo field_info;
  field_info.set_info("test", 0);
  field_info.set_type(std::type_index(typeid(int)), BCF_HT_INT);
  field_info.set_vcf_type(std::type_index(typeid(int)), BCF_HT_INT);
  auto& length_descriptor = field_info.m_length_descriptor;
  length_descriptor.resize(2u);
  //Length descriptors for each dimension
  length_descriptor.set_length_descriptor(0u, BCF_VL_VAR);
  length_descriptor.set_length_descriptor(1u, BCF_VL_VAR);
  //Delimiters
  length_descriptor.set_vcf_delimiter(0u, "|");
  length_descriptor.set_vcf_delimiter(1u, ",");
  GenomicsDBMultiDVectorField two_d_vector(field_info);
  std::string value = "1,2|3,4,5";
  two_d_vector.parse_and_store_numeric<int>(value.c_str(), value.length());
  //Index
  GenomicsDBMultiDVectorIdx dim_0_idx(&(two_d_vector.get_rw_data()[0]), &field_info);
  EXPECT_EQ(dim_0_idx.get_current_dim_index(), -1);
  //A[0]
  dim_0_idx.advance_to_index_in_next_dimension(0u);
  EXPECT_EQ(dim_0_idx.get_current_dim_index(), 0);
  EXPECT_EQ(dim_0_idx.get_current_index_in_current_dimension(), 0u);
  EXPECT_EQ(dim_0_idx.get_num_entries_in_current_dimension(), 2u);
  EXPECT_EQ(dim_0_idx.get_size_of_current_index(), 8u);
  auto data_ptr = dim_0_idx.get_ptr<int>();
  EXPECT_EQ(data_ptr[0u], 1);
  EXPECT_EQ(data_ptr[1u], 2);
  {
    //A[0][i]
    auto dim_1_idx = dim_0_idx;
    //A[0][0]
    dim_1_idx.advance_to_index_in_next_dimension(0u);
    EXPECT_EQ(dim_1_idx.get_current_dim_index(), 1);
    EXPECT_EQ(dim_1_idx.get_current_index_in_current_dimension(), 0u);
    EXPECT_EQ(dim_1_idx.get_num_entries_in_current_dimension(), 2u);
    EXPECT_EQ(dim_1_idx.get_size_of_current_index(), sizeof(int));
    EXPECT_EQ(dim_1_idx.get_element<int>(), 1);
    //A[0][1]
    dim_1_idx.advance_index_in_current_dimension();
    EXPECT_EQ(dim_1_idx.get_current_dim_index(), 1);
    EXPECT_EQ(dim_1_idx.get_current_index_in_current_dimension(), 1u);
    EXPECT_EQ(dim_1_idx.get_num_entries_in_current_dimension(), 2u);
    EXPECT_EQ(dim_1_idx.get_size_of_current_index(), sizeof(int));
    EXPECT_EQ(dim_1_idx.get_element<int>(), 2);
  }
  //A[1]
  dim_0_idx.advance_index_in_current_dimension();
  EXPECT_EQ(dim_0_idx.get_current_dim_index(), 0);
  EXPECT_EQ(dim_0_idx.get_current_index_in_current_dimension(), 1u);
  EXPECT_EQ(dim_0_idx.get_num_entries_in_current_dimension(), 2u);
  EXPECT_EQ(dim_0_idx.get_size_of_current_index(), 12u);
  data_ptr = dim_0_idx.get_ptr<int>();
  EXPECT_EQ(data_ptr[0u], 3);
  EXPECT_EQ(data_ptr[1u], 4);
  EXPECT_EQ(data_ptr[2u], 5);
  //VCF printer
  std::ostringstream fptr;
  GenomicsDBMultiDVectorFieldVCFPrinter print_op(fptr, field_info);
  two_d_vector.run_operation(print_op, &(two_d_vector.get_rw_data()[0]));
  EXPECT_EQ(fptr.str(), value);
}

TEST(multid_vector, 3D_test)
{
  FieldInfo field_info;
  field_info.set_info("test", 0);
  field_info.set_type(std::type_index(typeid(int)), BCF_HT_INT);
  field_info.set_vcf_type(std::type_index(typeid(int)), BCF_HT_INT);
  auto& length_descriptor = field_info.m_length_descriptor;
  length_descriptor.resize(3u);
  //Length descriptors for each dimension
  length_descriptor.set_length_descriptor(0u, BCF_VL_VAR);
  length_descriptor.set_length_descriptor(1u, BCF_VL_VAR);
  length_descriptor.set_length_descriptor(2u, BCF_VL_VAR);
  //Delimiters
  length_descriptor.set_vcf_delimiter(0u, "$");
  length_descriptor.set_vcf_delimiter(1u, "|");
  length_descriptor.set_vcf_delimiter(2u, ",");
  GenomicsDBMultiDVectorField three_d_vector(field_info);
  std::string value = "1,2|3,4,5$6,7,8,9|10,11|12";
  three_d_vector.parse_and_store_numeric<int>(value.c_str(), value.length());
  //Index
  GenomicsDBMultiDVectorIdx dim_0_idx(&(three_d_vector.get_rw_data()[0]), &field_info);
  EXPECT_EQ(dim_0_idx.get_current_dim_index(), -1);
  //A[0]
  dim_0_idx.advance_to_index_in_next_dimension(0u);
  EXPECT_EQ(dim_0_idx.get_current_dim_index(), 0);
  EXPECT_EQ(dim_0_idx.get_current_index_in_current_dimension(), 0u);
  EXPECT_EQ(dim_0_idx.get_num_entries_in_current_dimension(), 2u);
  EXPECT_EQ(dim_0_idx.get_size_of_current_index(), 60u);
  {
    //A[0][i]
    auto dim_1_idx = dim_0_idx;
    //A[0][0]
    dim_1_idx.advance_to_index_in_next_dimension(0u);
    EXPECT_EQ(dim_1_idx.get_current_dim_index(), 1);
    EXPECT_EQ(dim_1_idx.get_current_index_in_current_dimension(), 0u);
    EXPECT_EQ(dim_1_idx.get_num_entries_in_current_dimension(), 2u);
    EXPECT_EQ(dim_1_idx.get_size_of_current_index(), 8u);
    auto data_ptr = dim_1_idx.get_ptr<int>();
    EXPECT_EQ(data_ptr[0], 1);
    EXPECT_EQ(data_ptr[1], 2);
    //A[0][1]
    dim_1_idx.advance_index_in_current_dimension();
    EXPECT_EQ(dim_1_idx.get_current_dim_index(), 1);
    EXPECT_EQ(dim_1_idx.get_current_index_in_current_dimension(), 1u);
    EXPECT_EQ(dim_1_idx.get_num_entries_in_current_dimension(), 2u);
    EXPECT_EQ(dim_1_idx.get_size_of_current_index(), 12u);
    data_ptr = dim_1_idx.get_ptr<int>();
    EXPECT_EQ(data_ptr[0], 3);
    EXPECT_EQ(data_ptr[1], 4);
    EXPECT_EQ(data_ptr[2], 5);
    {
      //A[0][1][i]
      auto dim_2_idx = dim_1_idx;
      dim_2_idx.advance_to_index_in_next_dimension(0u);
      //A[0][1][0]
      EXPECT_EQ(dim_2_idx.get_current_dim_index(), 2);
      EXPECT_EQ(dim_2_idx.get_current_index_in_current_dimension(), 0u);
      EXPECT_EQ(dim_2_idx.get_num_entries_in_current_dimension(), 3u);
      EXPECT_EQ(dim_2_idx.get_size_of_current_index(), sizeof(int));
      EXPECT_EQ(dim_2_idx.get_element<int>(), 3);
      //A[0][1][1]
      dim_2_idx.advance_index_in_current_dimension();
      EXPECT_EQ(dim_2_idx.get_current_dim_index(), 2);
      EXPECT_EQ(dim_2_idx.get_current_index_in_current_dimension(), 1u);
      EXPECT_EQ(dim_2_idx.get_num_entries_in_current_dimension(), 3u);
      EXPECT_EQ(dim_2_idx.get_size_of_current_index(), sizeof(int));
      EXPECT_EQ(dim_2_idx.get_element<int>(), 4);
      //A[0][1][2]
      dim_2_idx.advance_index_in_current_dimension();
      EXPECT_EQ(dim_2_idx.get_current_dim_index(), 2);
      EXPECT_EQ(dim_2_idx.get_current_index_in_current_dimension(), 2u);
      EXPECT_EQ(dim_2_idx.get_num_entries_in_current_dimension(), 3u);
      EXPECT_EQ(dim_2_idx.get_size_of_current_index(), sizeof(int));
      EXPECT_EQ(dim_2_idx.get_element<int>(), 5);
    }
  }
  //A[1]
  dim_0_idx.advance_index_in_current_dimension();
  EXPECT_EQ(dim_0_idx.get_current_dim_index(), 0);
  EXPECT_EQ(dim_0_idx.get_current_index_in_current_dimension(), 1u);
  EXPECT_EQ(dim_0_idx.get_num_entries_in_current_dimension(), 2u);
  EXPECT_EQ(dim_0_idx.get_size_of_current_index(), 76u);
  {
    //A[1][i]
    auto dim_1_idx = dim_0_idx;
    //A[1][0]
    dim_1_idx.advance_to_index_in_next_dimension(0u);
    EXPECT_EQ(dim_1_idx.get_current_dim_index(), 1);
    EXPECT_EQ(dim_1_idx.get_current_index_in_current_dimension(), 0u);
    EXPECT_EQ(dim_1_idx.get_num_entries_in_current_dimension(), 3u);
    EXPECT_EQ(dim_1_idx.get_size_of_current_index(), 16u);
    auto data_ptr = dim_1_idx.get_ptr<int>();
    EXPECT_EQ(data_ptr[0], 6);
    EXPECT_EQ(data_ptr[1], 7);
    EXPECT_EQ(data_ptr[2], 8);
    EXPECT_EQ(data_ptr[3], 9);
    //A[0][1]
    dim_1_idx.advance_index_in_current_dimension();
    EXPECT_EQ(dim_1_idx.get_current_dim_index(), 1);
    EXPECT_EQ(dim_1_idx.get_current_index_in_current_dimension(), 1u);
    EXPECT_EQ(dim_1_idx.get_num_entries_in_current_dimension(), 3u);
    EXPECT_EQ(dim_1_idx.get_size_of_current_index(), 8u);
    data_ptr = dim_1_idx.get_ptr<int>();
    EXPECT_EQ(data_ptr[0], 10);
    EXPECT_EQ(data_ptr[1], 11);
    //A[0][2]
    dim_1_idx.advance_index_in_current_dimension();
    EXPECT_EQ(dim_1_idx.get_current_dim_index(), 1);
    EXPECT_EQ(dim_1_idx.get_current_index_in_current_dimension(), 2u);
    EXPECT_EQ(dim_1_idx.get_num_entries_in_current_dimension(), 3u);
    EXPECT_EQ(dim_1_idx.get_size_of_current_index(), 4u);
    data_ptr = dim_1_idx.get_ptr<int>();
    EXPECT_EQ(data_ptr[0], 12);
  }
  //VCF printer
  std::ostringstream fptr;
  GenomicsDBMultiDVectorFieldVCFPrinter print_op(fptr, field_info);
  three_d_vector.run_operation(print_op, &(three_d_vector.get_rw_data()[0]));
  EXPECT_EQ(fptr.str(), value);
}

TEST(multid_vector, 2D_test_with_missing_values)
{
  FieldInfo field_info;
  field_info.set_info("test", 0);
  field_info.set_type(std::type_index(typeid(int)), BCF_HT_INT);
  field_info.set_vcf_type(std::type_index(typeid(int)), BCF_HT_INT);
  auto& length_descriptor = field_info.m_length_descriptor;
  length_descriptor.resize(2u);
  //Length descriptors for each dimension
  length_descriptor.set_length_descriptor(0u, BCF_VL_VAR);
  length_descriptor.set_length_descriptor(1u, BCF_VL_VAR);
  //Delimiters
  length_descriptor.set_vcf_delimiter(0u, "|");
  length_descriptor.set_vcf_delimiter(1u, ",");
  GenomicsDBMultiDVectorField two_d_vector(field_info);
  std::string value = "|3,,NaN";
  two_d_vector.parse_and_store_numeric<int>(value.c_str(), value.length());
  //Index
  GenomicsDBMultiDVectorIdx dim_0_idx(&(two_d_vector.get_rw_data()[0]), &field_info);
  EXPECT_EQ(dim_0_idx.get_current_dim_index(), -1);
  //A[0]
  dim_0_idx.advance_to_index_in_next_dimension(0u);
  EXPECT_EQ(dim_0_idx.get_current_dim_index(), 0);
  EXPECT_EQ(dim_0_idx.get_current_index_in_current_dimension(), 0u);
  EXPECT_EQ(dim_0_idx.get_num_entries_in_current_dimension(), 2u);
  EXPECT_EQ(dim_0_idx.get_size_of_current_index(), 0u);
  //A[1]
  dim_0_idx.advance_index_in_current_dimension();
  EXPECT_EQ(dim_0_idx.get_current_dim_index(), 0);
  EXPECT_EQ(dim_0_idx.get_current_index_in_current_dimension(), 1u);
  EXPECT_EQ(dim_0_idx.get_num_entries_in_current_dimension(), 2u);
  EXPECT_EQ(dim_0_idx.get_size_of_current_index(), 3*sizeof(int));
  auto data_ptr = dim_0_idx.get_ptr<int>();
  EXPECT_EQ(data_ptr[0u], 3);
  EXPECT_EQ(data_ptr[1u], get_bcf_missing_value<int>());
  EXPECT_EQ(data_ptr[2u], get_bcf_missing_value<int>());
  //VCF printer
  std::ostringstream fptr;
  GenomicsDBMultiDVectorFieldVCFPrinter print_op(fptr, field_info);
  two_d_vector.run_operation(print_op, &(two_d_vector.get_rw_data()[0]));
  std::string serialization_return_value = "|3,,"; //NaN gets replaced with missing
  EXPECT_EQ(fptr.str(), serialization_return_value);
}
