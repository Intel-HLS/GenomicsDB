#ifdef NDEBUG	//needs asserts
#undef NDEBUG
#endif

#include <iostream>
#include<fstream>
#include<unordered_map>
#include <vector>
#include <sstream>
#include "lut.h"

#include "vcf.h"

using namespace std;

int main(int argc, char** argv) {
  ifstream gold;
  ifstream test;
  if(argc < 3)
    exit(-1);
  gold.open(argv[1]);
  test.open(argv[2]);
  assert(gold.is_open() && test.is_open());

  string gold_line;
  string test_line;
  CombineAllelesLUT alleles_LUT(2);
  vector<string> gold_tokens;
  vector<string> test_tokens;
  uint64_t line_num = 0ull;
  while(1)
  {
    istringstream iss_gold;
    istringstream iss_test;
    
    gold_tokens.resize(100);
    test_tokens.resize(100);


    unsigned idx = 0u;
    unsigned num_alleles = 0u;

    getline(gold, gold_line);
    getline(test, test_line);
    //same number of lines 
    assert((gold.eof() && test.eof())  || (!gold.eof() && !test.eof()));
    if(gold.eof())
      break;
    
    iss_gold.str(gold_line);
    while(!iss_gold.eof())
    {
      getline(iss_gold, gold_tokens[idx],',');
      if(gold_tokens[idx] == "<NON_REF>")
        num_alleles = idx - 1;
      ++idx;
      if(idx >= gold_tokens.size())
	gold_tokens.resize(2*idx+1);
    }
    //Erase first 2 elements
    gold_tokens.erase(gold_tokens.begin());
    gold_tokens.erase(gold_tokens.begin());

    iss_test.str(test_line);
    idx = 0u;
    while(!iss_test.eof())
    {
      getline(iss_test, test_tokens[idx],',');
      if(test_tokens[idx] == "<NON_REF>")
        assert(num_alleles == idx);     //Same #alleles
      ++idx;
      if(idx >= test_tokens.size())
	test_tokens.resize(2*idx+1);
    }
    int64_t position = strtoll(test_tokens[0].c_str(), 0, 10);
    //Erase first element
    test_tokens.erase(test_tokens.begin());
    //Same REF
    assert(test_tokens[0].length() == gold_tokens[0].length());
    if(test_tokens[0] != gold_tokens[0])
    {
      assert(test_tokens[0].length() == 1 && test_tokens[0][0] == 'N');
      test_tokens[0] = gold_tokens[0];
    }
    alleles_LUT.reset_luts();
    alleles_LUT.resize_luts_if_needed(num_alleles);
    unordered_map<string, unsigned> gold_allele_2_idx;
    for(auto i=0u;i<num_alleles;++i)
      gold_allele_2_idx[gold_tokens[i]] = i;
    for(auto i=0u;i<num_alleles;++i)
    {
      auto iter = gold_allele_2_idx.find(test_tokens[i]);
      assert(iter != gold_allele_2_idx.end());
      alleles_LUT.add_input_merged_idx_pair(1u, i, (*iter).second);
    }
    unsigned num_gts = (num_alleles*(num_alleles+1))/2;
    for(auto i=0u;i<num_alleles;++i)
    {
      auto test_i = alleles_LUT.get_input_idx_for_merged(1u, i);
      for(auto j=i;j<num_alleles;++j)
      {
        auto test_j = alleles_LUT.get_input_idx_for_merged(1u, j);
        assert(test_tokens[num_alleles + bcf_alleles2gt(test_i, test_j)] == gold_tokens[num_alleles + bcf_alleles2gt(i, j)]);
      }
    }
    ++line_num;
  }
  gold.close();
  test.close();
  return 0;
}
