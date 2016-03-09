#ifndef COMMON_HEADERS_H
#define COMMON_HEADERS_H

#include <stdio.h>
#include <typeinfo>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sys/stat.h>
#include <assert.h>
#include <math.h>
#include <queue>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <set>
#include <map>
#include <climits>
#include <sstream>
#include <iomanip>
#include <exception>
#include <fstream>
#include <functional>
#include <omp.h>
#include <typeindex>

typedef std::pair<int64_t, int64_t> ColumnRange;
typedef std::pair<int64_t, int64_t> RowRange;
bool ColumnRangeCompare(const ColumnRange& x, const ColumnRange& y);

class CircularBufferController
{
  public:
    CircularBufferController(unsigned num_entries)
    {
      m_num_entries = num_entries;
      m_num_entries_with_valid_data = 0u;
      m_num_reserved_entries = 0u;
      m_curr_write_idx = 0u;
      m_curr_read_idx = 0u;
    }
    //Advance write idx - increase #valid entries, decrease #reserved entries if specified
    inline void advance_write_idx(bool unreserve=false)
    {
      m_curr_write_idx = (m_curr_write_idx+1u)%m_num_entries;
      assert(m_num_entries_with_valid_data < m_num_entries);
      ++m_num_entries_with_valid_data;
      if(unreserve)
      {
        assert(m_num_reserved_entries > 0u);
        --m_num_reserved_entries;
      }
    }
    //Reserves an entry without marking it as valid
    void reserve_entry()   { ++m_num_reserved_entries; }
    inline void advance_read_idx()
    {
      m_curr_read_idx = (m_curr_read_idx+1u)%m_num_entries;
      assert(m_num_entries_with_valid_data > 0u);
      --m_num_entries_with_valid_data;
    }
    inline unsigned get_num_entries_with_valid_data() const
    { return m_num_entries_with_valid_data; }
    inline unsigned get_num_empty_entries() const
    { return m_num_entries - m_num_entries_with_valid_data - m_num_reserved_entries; }
    //Get idx
    inline unsigned get_write_idx() const { return m_curr_write_idx; }
    inline unsigned get_read_idx() const { return m_curr_read_idx; }
  protected:
    //Points to latest entry with valid data
    unsigned m_curr_write_idx;
    //Points to entry being read currently
    unsigned m_curr_read_idx;
    unsigned m_num_entries;
    unsigned m_num_entries_with_valid_data;
    unsigned m_num_reserved_entries;
};


#endif
