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

#include <test_mapping_data_loader.h>

MappingDataLoaderTester::~MappingDataLoaderTester() {
  m_contig_idx_to_info.clear();
  m_contig_begin_2_idx.clear();
  m_contig_end_2_idx.clear();
  m_row_idx_to_info.clear();
  m_callset_name_to_row_idx.clear();
  m_field_name_to_idx.clear();
  m_field_idx_to_info.clear();
}

MappingDataLoaderTester::MappingDataLoaderTester(const SQLVidMapperRequest& request) {
  SQLBasedVidMapper sql_mapper(request);
  sql_mapper.load_mapping_data_from_db();
  m_contig_idx_to_info = sql_mapper.get_contigs();
  m_contig_begin_2_idx = sql_mapper.get_contig_begin_offsets();
  m_contig_end_2_idx = sql_mapper.get_contig_end_offsets();
  m_row_idx_to_info = sql_mapper.get_callsets();
  m_callset_name_to_row_idx = sql_mapper.get_callset_name_to_idx_map();
  m_field_name_to_idx = sql_mapper.get_field_name_to_idx_map();
  m_field_idx_to_info = sql_mapper.get_fields();
}

void MappingDataLoaderTester::validate_contig_info() {
  std::cout <<"------------------------------------------------------\n";
  int64_t prev_offset = -1;

  for (auto& contig : m_contig_idx_to_info) {
    std::cout <<contig.m_contig_idx <<" - " <<contig.m_name <<" - " 
      <<contig.m_tiledb_column_offset <<" - " <<contig.m_length <<"\n";
    EXPECT_GT(contig.m_tiledb_column_offset, prev_offset);
    prev_offset = contig.m_tiledb_column_offset;
  }

  std::cout <<"------------------------------------------------------\n";
  int64_t prev_begin_first = -1;

  for (auto& contig_begin : m_contig_begin_2_idx) {
    std::cout <<contig_begin.first <<" - " <<contig_begin.second <<"\n";
    EXPECT_GT(contig_begin.first, prev_begin_first);
    prev_begin_first = contig_begin.first;
  }

  std::cout <<"------------------------------------------------------\n";
  int64_t prev_end_first = -1;

  for (auto& contig_end : m_contig_end_2_idx) {
    std::cout <<contig_end.first <<" - " <<contig_end.second <<"\n";
    EXPECT_GT(contig_end.first, prev_end_first);
    prev_end_first = contig_end.first;
  }

  std::cout <<"------------------------------------------------------\n";
  return;
}

void MappingDataLoaderTester::validate_callset_info() {
  std::cout <<"------------------------------------------------------\n";
  int64_t index = -1;

  for (std::vector<CallSetInfo>::iterator it = m_row_idx_to_info.begin(); it != m_row_idx_to_info.end(); ++it) {
    index++;
    std::cout <<it->m_row_idx <<" - " <<it->m_name <<" - " <<it->m_file_idx <<" - " <<it->m_idx_in_file <<"\n";
    EXPECT_EQ(index, it->m_row_idx);
    EXPECT_EQ(m_callset_name_to_row_idx[it->m_name], index);
  }

  std::cout <<"------------------------------------------------------\n";
  return;
}

void MappingDataLoaderTester::validate_field_info() {
  std::cout <<"------------------------------------------------------\n";
  int64_t index = -1;

  for (auto& field : m_field_idx_to_info) {
    std::cout <<field.m_field_idx <<" - " <<field.m_name <<" - " 
      <<field.m_length_descriptor <<" - " <<field.m_num_elements <<" - " 
      <<field.m_bcf_ht_type <<" - " <<field.m_is_vcf_FILTER_field <<" - " 
      <<field.m_is_vcf_INFO_field <<" - " <<field.m_is_vcf_FORMAT_field <<" - "
      <<field.m_VCF_field_combine_operation <<"\n";
    EXPECT_EQ(++index, field.m_field_idx);
    EXPECT_EQ(m_field_name_to_idx[field.m_name], index);
  }

  std::cout <<"------------------------------------------------------\n";
  return;
}

MappingDataLoaderTester* SQLMapperTest::loaderTester = NULL;

TEST_F(SQLMapperTest, ValidContigInfo) {
  loaderTester->validate_contig_info();
}

TEST_F(SQLMapperTest, ValidCallSetInfo) {
  loaderTester->validate_callset_info();
}

TEST_F(SQLMapperTest, ValidFieldInfo) {
  loaderTester->validate_field_info();
}
