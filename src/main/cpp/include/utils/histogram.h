/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

#ifndef HISTOGRAM_H
#define  HISTOGRAM_H

#include "headers.h"

class HistogramException : public std::exception {
  public:
    HistogramException(const std::string m="") : msg_("HistogramException : "+m) { ; }
    ~HistogramException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class UnimplementedHistogramOperation : public std::exception {
  public:
    UnimplementedHistogramOperation(const std::string m="") : msg_("Unimplemented histogram operation"+m) { ; }
    ~UnimplementedHistogramOperation() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class Histogram
{
  public:
    Histogram()
    {
      m_total = 0;
      m_lo = 0;
      m_hi = 0;
      clear();
    }
    ~Histogram() { clear(); }
    void clear();
    class HistogramIterator
    {
      public:
        HistogramIterator(const Histogram* ptr, unsigned bin_idx, unsigned num_bins)
        {
          m_histogram_ptr = ptr;
          m_bin_idx = bin_idx;
          m_num_bins = num_bins;
        }
        bool operator!=(const HistogramIterator& other)
        {
          return m_histogram_ptr != other.m_histogram_ptr || m_bin_idx != other.m_bin_idx;
        }
        const HistogramIterator& operator++()
        {
          ++m_bin_idx;
          return *this;
        }
        uint64_t operator*() const { return m_histogram_ptr->get_histogram_value(m_bin_idx); }
        uint64_t get_lo() const { return m_histogram_ptr->get_lo(m_bin_idx); }
        uint64_t get_hi() const { return m_histogram_ptr->get_hi(m_bin_idx); }
        unsigned get_bin_idx() const { return m_bin_idx; }
      private:
        unsigned m_bin_idx;
        unsigned m_num_bins;
        const Histogram* m_histogram_ptr;
    };
    typedef const HistogramIterator const_iterator;
    typedef HistogramIterator iterator;
    uint64_t get_histogram_value(unsigned bin_idx) const
    {
      assert(bin_idx < m_histogram_bins.size());
      return m_histogram_bins[bin_idx];
    }
    iterator begin() { return HistogramIterator(this, 0u, m_histogram_bins.size()); }
    iterator end() { return HistogramIterator(this, m_histogram_bins.size(), m_histogram_bins.size()); }
    const_iterator begin() const { return HistogramIterator(this, 0u, m_histogram_bins.size()); }
    const_iterator end() const { return HistogramIterator(this, m_histogram_bins.size(), m_histogram_bins.size()); }
    uint64_t get_total() const { return m_total; }
    void add_value(uint64_t value);
    void add_interval(uint64_t lo, uint64_t hi);
    void print(std::ostream& fptr=std::cout) const;
    void reset_counters();
    //For MPI communication
    uint64_t serialize(uint8_t*& data, uint64_t curr_offset, bool realloc_if_needed=true) const;
    uint64_t deserialize(const uint8_t* data, uint64_t offset);
    //Abstract functions - must override 
    virtual unsigned get_bin_idx_for_value(uint64_t value) const = 0; 
    virtual uint64_t get_lo(unsigned bin_idx) const = 0;
    virtual uint64_t get_hi(unsigned bin_idx) const = 0;
  protected:
    std::vector<uint64_t> m_histogram_bins;
    uint64_t m_total;
    uint64_t m_lo;
    uint64_t m_hi;
};

class UniformHistogram : public Histogram
{
  public:
    //[lo, hi)
    UniformHistogram(uint64_t lo, uint64_t hi, unsigned num_intervals);
    virtual ~UniformHistogram() = default;
    unsigned get_bin_idx_for_value(uint64_t value) const 
    {
      unsigned bin_idx = (value - m_lo)/m_size_of_bin;
      assert(bin_idx < m_histogram_bins.size());
      return bin_idx; 
    }
    uint64_t get_lo(unsigned bin_idx) const
    {
      assert(bin_idx < m_histogram_bins.size());
      return m_lo + bin_idx*m_size_of_bin;
    }
    uint64_t get_hi(unsigned bin_idx) const
    {
      assert(bin_idx < m_histogram_bins.size());
      return m_lo + (bin_idx+1)*m_size_of_bin - 1;
    }
    void sum_up_histogram(const UniformHistogram& other);
    void sum_up_histogram(const UniformHistogram* other) { sum_up_histogram(*other); }
    //For MPI communication
    uint64_t serialize(uint8_t*& data, uint64_t curr_offset, bool realloc_if_needed=true) const;
    uint64_t deserialize(const uint8_t* data, uint64_t offset);
  private:
    uint64_t m_size_of_bin;
};

#endif
