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
}

MappingDataLoaderTester::MappingDataLoaderTester(const SQLVidMapperRequest& request) {
  SQLBasedVidMapper sql_mapper(request);
  sql_mapper.load_mapping_data_from_db();
  m_contig_idx_to_info = sql_mapper.get_contigs();
  m_contig_begin_2_idx = sql_mapper.get_contig_begin_offsets();
  m_contig_end_2_idx = sql_mapper.get_contig_end_offsets();
  m_row_idx_to_info = sql_mapper.get_callsets();
  m_callset_name_to_row_idx = sql_mapper.get_name_to_idx_map();
}

void MappingDataLoaderTester::validate_contig_info() {
  std::cout <<"------------------------------------------------------\n";
  int64_t prev_offset = -1;

  for (std::vector<ContigInfo>::iterator it = m_contig_idx_to_info.begin(); it != m_contig_idx_to_info.end(); ++it) {
    std::cout <<it->m_contig_idx <<" - " <<it->m_name <<" - " <<it->m_tiledb_column_offset <<" - " <<it->m_length <<"\n";
    EXPECT_GT(it->m_tiledb_column_offset, prev_offset);
    prev_offset = it->m_tiledb_column_offset;
  }

  std::cout <<"------------------------------------------------------\n";
  int64_t prev_begin_first = -1;

  for (std::vector<std::pair<int64_t, int>>::iterator it = m_contig_begin_2_idx.begin(); it != m_contig_begin_2_idx.end(); ++it) {
    std::cout <<it->first <<" - " <<it->second <<"\n";
    EXPECT_GT(it->first, prev_begin_first);
    prev_begin_first = it->first;
  }

  std::cout <<"------------------------------------------------------\n";
  int64_t prev_end_first = -1;

  for (std::vector<std::pair<int64_t, int>>::iterator it = m_contig_end_2_idx.begin(); it != m_contig_end_2_idx.end(); ++it) {
    std::cout <<it->first <<" - " <<it->second <<"\n";
    EXPECT_GT(it->first, prev_end_first);
    prev_end_first = it->first;
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

MappingDataLoaderTester* SQLMapperTest::loaderTester = NULL;

TEST_F(SQLMapperTest, ValidContigInfo) {
  loaderTester->validate_contig_info();
}

TEST_F(SQLMapperTest, ValidCallSetInfo) {
  loaderTester->validate_callset_info();
}
