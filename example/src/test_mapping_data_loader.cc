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
#include <getopt.h>
#include "headers.h"
#include "vid_mapper_sql.h"
#include <cassert>

const std::string HOSTNAME_PARAM = "host_name";
const std::string HOSTNAME_DEFAULT = "localhost";
const std::string USERNAME_PARAM = "user_name";
const std::string USERNAME_DEFAULT = "postgres";
const std::string PASSWD_PARAM = "pass_word";
const std::string PASSWD_DEFAULT = "postgres";
const std::string DBNAME_PARAM = "db_name";
const std::string DBNAME_DEFAULT = "gendb";
const std::string WORKSPACE_PARAM = "work_space";
const std::string WORKSPACE_DEFAULT = "/home/rmantrix/git_repos/NomixDB/INST_GenomicsDB/tests/workspace";
const std::string ARRAYNAME_PARAM = "array_name";
const std::string ARRAYNAME_DEFAULT = "rsm_test";

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
  SQLBasedVidMapper vid_mapper(request);
  vid_mapper.load_mapping_data_from_db();
  m_contig_idx_to_info = vid_mapper.get_contigs();
  m_contig_begin_2_idx = vid_mapper.get_contig_begin_offsets();
  m_contig_end_2_idx = vid_mapper.get_contig_end_offsets();
  m_row_idx_to_info = vid_mapper.get_callsets();
  m_callset_name_to_row_idx = vid_mapper.get_name_to_idx_map();
}

void MappingDataLoaderTester::validate_contig_info() {
  std::cout <<"------------------------------------------------------\n";
  int64_t prev_offset = -1;

  for (std::vector<ContigInfo>::iterator it = m_contig_idx_to_info.begin(); it != m_contig_idx_to_info.end(); ++it) {
    std::cout <<it->m_contig_idx <<" - " <<it->m_name <<" - " <<it->m_tiledb_column_offset <<" - " <<it->m_length <<"\n";
    assert(it->m_tiledb_column_offset > prev_offset);
    prev_offset = it->m_tiledb_column_offset;
  }

  std::cout <<"------------------------------------------------------\n";
  int64_t prev_begin_first = -1;

  for (std::vector<std::pair<int64_t, int>>::iterator it = m_contig_begin_2_idx.begin(); it != m_contig_begin_2_idx.end(); ++it) {
    std::cout <<it->first <<" - " <<it->second <<"\n";
    assert(it->first > prev_begin_first);
    prev_begin_first = it->first;
  }

  std::cout <<"------------------------------------------------------\n";
  int64_t prev_end_first = -1;

  for (std::vector<std::pair<int64_t, int>>::iterator it = m_contig_end_2_idx.begin(); it != m_contig_end_2_idx.end(); ++it) {
    std::cout <<it->first <<" - " <<it->second <<"\n";
    assert(it->first > prev_end_first);
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
    assert(index == it->m_row_idx);
    assert(m_callset_name_to_row_idx[it->m_name] == index);
  }

  std::cout <<"------------------------------------------------------\n";
  return;
}

void get_input(const std::string param_name, const std::string param_value, std::string& location) {
  std::string value;
  std::cout <<"Enter value for <" <<param_name <<"> - type DEFAULT to use default value: ";
  std::cin >> value;
  if (0 == value.compare("DEFAULT")) {
      location = param_value;
  } else {
      location = value;
  }
  return;
}

int main(int argc, char *argv[]) {
  std::cout <<"------------ MappingData Loader - BEGIN -----------------\n";
  SQLVidMapperRequest request;

  get_input(HOSTNAME_PARAM, HOSTNAME_DEFAULT, request.host_name);
  get_input(USERNAME_PARAM, postgres, request.user_name);
  get_input(PASSWD_PARAM, PASSWD_DEFAULT, request.pass_word);
  get_input(DBNAME_PARAM, DBNAME_DEFAULT, request.db_name);

  get_input(WORKSPACE_PARAM, WORKSPACE_DEFAULT, request.work_space);
  //request.array_name = "hg19";
  get_input(ARRAYNAME_PARAM, ARRAYNAME_DEFAULT, request.array_name);

  MappingDataLoaderTester tester(request);

  tester.validate_contig_info();

  tester.validate_callset_info();

  std::cout <<"------------ MappingData Loader -  END  -----------------\n";
  return(0);
}

