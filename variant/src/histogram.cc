#include "histogram.h"

void Histogram::clear()
{
  m_histogram_bins.clear();
}

void Histogram::print(std::ostream& fptr) const
{
  fptr << "Histogram: { "; 
  for(auto b=this->begin(), e = this->end();b!=e;++b)
  {
    if((*b) != 0)
      fptr << "[" << b.get_lo() <<","<<b.get_hi()<<"]:"<<*b << ", ";
  }
  fptr << " }\n";
}

void Histogram::add_value(uint64_t value)
{
  if(value < m_lo || value > m_hi)
    throw HistogramException("Value "+std::to_string(value)+" is outsize the range of histogram ["+std::to_string(m_lo)
        +","+std::to_string(m_hi)+"]"); 
  unsigned bin_idx = get_bin_idx_for_value(value);
  ++(m_histogram_bins[bin_idx]);
  ++m_total;
}

void Histogram::add_interval(uint64_t lo, uint64_t hi)
{
  if(lo == hi)
  {
    add_value(lo);
    return;
  }
  if(lo < m_lo || lo > m_hi)
    throw HistogramException("Value "+std::to_string(lo)+" is outsize the range of histogram ["+std::to_string(m_lo)
        +","+std::to_string(m_hi)+"]"); 
  if(hi < m_lo || hi > m_hi)
    throw HistogramException("Value "+std::to_string(hi)+" is outsize the range of histogram ["+std::to_string(m_lo)
        +","+std::to_string(m_hi)+"]"); 
  unsigned lo_bin_idx = get_bin_idx_for_value(lo);
  unsigned hi_bin_idx = get_bin_idx_for_value(hi);
  for(auto i=lo_bin_idx;i<=hi_bin_idx;++i)
  {
    ++(m_histogram_bins[i]);
    ++m_total;
  }
}

uint64_t Histogram::serialize(uint8_t*& data, uint64_t curr_offset, bool realloc_if_needed) const
{
  uint64_t vec_size = m_histogram_bins.size()*sizeof(uint64_t);
  uint64_t to_add = 4*sizeof(uint64_t) + vec_size;
  auto new_size = curr_offset + to_add;
  if(realloc_if_needed)
    data = (uint8_t*)realloc(data, new_size);
  auto idx = curr_offset;
  uint64_t num_bins = m_histogram_bins.size();
  memcpy(&(data[idx]), &num_bins, sizeof(uint64_t));
  idx += sizeof(uint64_t);
  memcpy(&(data[idx]), &m_total, sizeof(uint64_t));
  idx += sizeof(uint64_t);
  memcpy(&(data[idx]), &m_lo, sizeof(uint64_t));
  idx += sizeof(uint64_t);
  memcpy(&(data[idx]), &m_hi, sizeof(uint64_t));
  idx += sizeof(uint64_t);
  memcpy(&(data[idx]), &(m_histogram_bins[0]), vec_size);
  idx += vec_size;
  return idx;
}

uint64_t Histogram::deserialize(const uint8_t* data, uint64_t offset)
{
  auto idx = offset;
  uint64_t num_bins = 0;
  memcpy(&num_bins, &(data[idx]), sizeof(uint64_t));
  idx += sizeof(uint64_t);
  memcpy(&m_total, &(data[idx]),sizeof(uint64_t));
  idx += sizeof(uint64_t);
  memcpy(&m_lo, &(data[idx]),sizeof(uint64_t));
  idx += sizeof(uint64_t);
  memcpy(&m_hi, &(data[idx]),sizeof(uint64_t));
  idx += sizeof(uint64_t);
  m_histogram_bins.resize(num_bins);
  uint64_t vec_size = m_histogram_bins.size()*sizeof(uint64_t);
  memcpy(&(m_histogram_bins[0]), &(data[idx]), vec_size);
  idx += vec_size;
  return idx;
}

void Histogram::reset_counters()
{
  for(auto i=0u;i<m_histogram_bins.size();++i)
    m_histogram_bins[i] = 0ull;
  m_total = 0ull;
}

UniformHistogram::UniformHistogram(uint64_t lo, uint64_t hi, unsigned num_intervals)
  : Histogram()
{
  m_lo = lo;
  m_hi = hi;
  m_size_of_bin = (hi - lo)/num_intervals;
  m_histogram_bins.resize(num_intervals);
  for(auto i=0u;i<num_intervals;++i)
    m_histogram_bins[i] = 0;
}

void UniformHistogram::sum_up_histogram(const UniformHistogram& other)
{
  Histogram::iterator curr = begin();
  for(auto b=other.begin(),e=other.end();b!=e;++b)
  {
    if(curr.get_lo() != b.get_lo() || curr.get_hi() != b.get_hi())
      throw HistogramException("To sum up UniformHistogram objects, bin ranges must match");
    auto value = *b;
    m_histogram_bins[curr.get_bin_idx()] += value;
    ++curr;
  }
  m_total += other.get_total();
}

uint64_t UniformHistogram::serialize(uint8_t*& data, uint64_t curr_offset, bool realloc_if_needed) const
{
  auto idx = Histogram::serialize(data, curr_offset, realloc_if_needed);
  auto new_size = idx + sizeof(uint64_t);
  if(realloc_if_needed)
    data = (uint8_t*)realloc(data, new_size);
  memcpy(&(data[idx]), &m_size_of_bin, sizeof(uint64_t));
  idx += sizeof(uint64_t);
  return idx;
}

uint64_t UniformHistogram::deserialize(const uint8_t* data, uint64_t offset)
{
  auto idx = Histogram::deserialize(data, offset);
  memcpy(&m_size_of_bin, &(data[idx]), sizeof(uint64_t));
  idx += sizeof(uint64_t);
  return idx;
}
