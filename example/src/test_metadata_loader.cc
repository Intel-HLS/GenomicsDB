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


class MetaDataLoaderTester {
  private:
    SQLBasedVidMapper vid_mapper;
    std::vector<ContigInfo> m_contig_idx_to_info;
    std::vector<std::pair<int64_t, int>> m_contig_begin_2_idx;
    std::vector<std::pair<int64_t, int>> m_contig_end_2_idx;
  public:
    MetaDataLoaderTester();
    void print_contig_info();
    ~MetaDataLoaderTester();
};

MetaDataLoaderTester::~MetaDataLoaderTester() {
  m_contig_idx_to_info.clear();
  m_contig_begin_2_idx.clear();
  m_contig_end_2_idx.clear();
}

MetaDataLoaderTester::MetaDataLoaderTester() {
  vid_mapper.load_metadata_from_db();
  m_contig_idx_to_info = vid_mapper.get_contigs();
  m_contig_begin_2_idx = vid_mapper.get_contig_begin();
  m_contig_end_2_idx = vid_mapper.get_contig_end();
}

void MetaDataLoaderTester::print_contig_info() {
  std::cout <<"------------------------------------------------------\n";

  for (std::vector<ContigInfo>::iterator it = m_contig_idx_to_info.begin(); it != m_contig_idx_to_info.end(); ++it) {
    std::cout <<it->m_contig_idx <<" - " <<it->m_name <<" - " <<it->m_tiledb_column_offset <<" - " <<it->m_length <<"\n";
  }

  std::cout <<"------------------------------------------------------\n";

  for (std::vector<std::pair<int64_t, int>>::iterator it = m_contig_begin_2_idx.begin(); it != m_contig_begin_2_idx.end(); ++it) {
    std::cout <<it->first <<" - " <<it->second <<"\n";
  }

  std::cout <<"------------------------------------------------------\n";

  for (std::vector<std::pair<int64_t, int>>::iterator it = m_contig_end_2_idx.end(); it != m_contig_end_2_idx.end(); ++it) {
    std::cout <<it->first <<" - " <<it->second <<"\n";
  }

  std::cout <<"------------------------------------------------------\n";
  return;
}

int main(int argc, char *argv[]) {
  std::cout <<"------------ MetaData Loader - BEGIN -----------------\n";

  MetaDataLoaderTester tester;

  tester.print_contig_info();

  std::cout <<"------------ MetaData Loader -  END  -----------------\n";
  return(0);
}
