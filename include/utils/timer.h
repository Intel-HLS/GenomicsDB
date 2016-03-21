/**
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

#ifndef TIMER_H
#define TIMER_H

#include<vector>
#include<stdlib.h>
#include<time.h>
#include<sys/time.h>
#include <iomanip>

class Timer
{
  public:
    Timer()
    { 
      m_cumulative_cpu_time = 0;
      m_cumulative_wall_clock_time = 0;
      start();
    };
    inline void start()
    {
      /*m_begin_cpu_time = clock();*/
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &m_begin_cpu_time);
      gettimeofday(&m_begin_wall_clock_time, 0);
    }
    inline void stop()
    {
      /*CPU time*/
      /*m_last_interval_cpu_time = ((double)(clock() - m_begin_cpu_time))/CLOCKS_PER_SEC;*/
      struct timespec end_cpu_time;
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end_cpu_time);
      /*Wall clock time*/
      struct timeval end_time;
      gettimeofday(&end_time, 0);
      m_last_interval_cpu_time = ((double)(static_cast<uint64_t>(end_cpu_time.tv_sec - m_begin_cpu_time.tv_sec)*1000000000ull
            + static_cast<uint64_t>(end_cpu_time.tv_nsec - m_begin_cpu_time.tv_nsec)))/1000000000ull;
      m_last_interval_wall_clock_time = ((double)((end_time.tv_sec - m_begin_wall_clock_time.tv_sec)*1000000ull
            + (end_time.tv_usec - m_begin_wall_clock_time.tv_usec)))/1000000;
      m_cumulative_cpu_time += m_last_interval_cpu_time;
      m_cumulative_wall_clock_time += m_last_interval_wall_clock_time;
    }
    void print(const std::string& prefix="", std::ostream& fptr = std::cout) const
    {
      if(prefix.size() > 0u)
        fptr << prefix <<" : ";
      fptr << "Wall-clock time(s) : "<< std::setprecision(6) << m_last_interval_wall_clock_time << " Cpu time(s) : "
        << m_last_interval_cpu_time << "\n";
    }
    void print_cumulative(const std::string& prefix="", std::ostream& fptr = std::cout) const
    {
      if(prefix.size() > 0u)
        fptr << prefix <<" : ";
      fptr << "Wall-clock time(s) : "<< std::setprecision(6) << m_cumulative_wall_clock_time << " Cpu time(s) : "
        << m_cumulative_cpu_time << "\n";
    }
    void get_last_interval_times(std::vector<double>& timings, unsigned timer_idx) const
    {
      auto idx = 2u*timer_idx;
      assert(idx + 2u <= timings.size());
      timings[idx] = m_last_interval_cpu_time;
      timings[idx+1] = m_last_interval_wall_clock_time;
    }
  private:
    /*clock_t m_begin_cpu_time;*/
    struct timespec m_begin_cpu_time;
    struct timeval m_begin_wall_clock_time;
    //seconds
    double m_last_interval_cpu_time;
    double m_last_interval_wall_clock_time;
    double m_cumulative_cpu_time;
    double m_cumulative_wall_clock_time;
};

#endif
