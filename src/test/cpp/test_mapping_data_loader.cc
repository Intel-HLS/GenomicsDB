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

#include <iostream>
#include <fstream>
#include <getopt.h>
#include "headers.h"
#include "vid_mapper_sql.h"
#include <limits.h>
#include "gtest/gtest.h"

const std::string config_file_name = "/tmp/sql_mapper/sql_mapper_config.txt";

class MappingDataLoaderTester {
  private:
    std::vector<ContigInfo> m_contig_idx_to_info;
    std::vector<std::pair<int64_t, int>> m_contig_begin_2_idx;
    std::vector<std::pair<int64_t, int>> m_contig_end_2_idx;
    std::vector<CallSetInfo> m_row_idx_to_info;
    std::unordered_map<std::string, int64_t> m_callset_name_to_row_idx;
  public:
    MappingDataLoaderTester(const SQLVidMapperRequest&);
    void validate_contig_info();
    void validate_callset_info();
    ~MappingDataLoaderTester();
};

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

/**
 * Following config parameters are expected in the file:
 * /tmp/sql_mapper/sql_mapper_config.txt
 * If the file does not exist, then default values are used
 * "host_name"  - "localhost"
 * "user_name" - "postgres"
 * "pass_word" - "postgres"
 * "db_name" - "gendb"
 * "work_space" - "/tmp/sql_mapper/workspace"
 * "array_name" - "sql_mapper_test"
 */
class SQLMapperTest : public ::testing::Test {
  public:
    static MappingDataLoaderTester* loaderTester;

  protected:
    static void SetUpTestCase() {
      std::cout <<"------------ MappingData Loader - BEGIN -----------------\n";
      SQLVidMapperRequest request;
      std::ifstream config_file(config_file_name);

      if (config_file.is_open()) {
        getline(config_file, request.host_name);
        getline(config_file, request.user_name);
        getline(config_file, request.pass_word);
        getline(config_file, request.db_name);
        getline(config_file, request.work_space);
        getline(config_file, request.array_name);
        config_file.close();
      } else {
        request.host_name = "localhost";
        request.user_name = "postgres";
        request.pass_word = "postgres";
        request.db_name = "gendb";
        request.work_space = "/tmp/sql_mapper/workspace";
        request.array_name = "sql_mapper_test";
      }

      loaderTester = new MappingDataLoaderTester(request);
      std::cout <<"------------ MappingData Loader -  END  -----------------\n";
    }

    static void TearDownTestCase() {
        delete loaderTester;
        loaderTester = NULL;
    }

    virtual void SetUp() {
    }

    virtual void TearDown() {
    }
};

MappingDataLoaderTester* SQLMapperTest::loaderTester = NULL;

TEST_F(SQLMapperTest, ValidContigInfo) {
  loaderTester->validate_contig_info();
}

TEST_F(SQLMapperTest, ValidCallSetInfo) {
  loaderTester->validate_callset_info();
}
