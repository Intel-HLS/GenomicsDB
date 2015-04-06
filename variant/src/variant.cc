#include "variant.h"

std::string g_non_reference_allele = "<NON_REF>";

GTColumn::GTColumn(int64_t col, uint64_t row_num) {
  ALT_.resize(row_num);
  col_ = col;
  REF_.resize(row_num);
  PL_.resize(row_num);
  AF_.resize(row_num);
  AN_.resize(row_num);
  AC_.resize(row_num);
}

void GTColumn::reset()
{
  for(auto i=0ull;i<REF_.size();++i)
  {
    REF_[i].clear();
    ALT_[i].clear();
    PL_[i].clear();
    AF_.clear();
    AN_.clear();
    AC_.clear();
  }
}
