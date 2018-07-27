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

#ifndef TIMER_H
#define TIMER_H

#include<vector>
#include<stdlib.h>
#include<time.h>
#include<sys/time.h>
#include <iomanip>

//MacOS does not have a clock_gettime function
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#define CLOCK_THREAD_CPUTIME_ID 0
int clock_gettime(clock_id_t clk_id, struct timespec *tp);
#endif

class Timer
{
  public:
    Timer()
    {
      m_last_interval_wall_clock_time = 0;
      m_last_interval_cpu_time = 0;
      m_cumulative_wall_clock_time = 0;
      m_cumulative_cpu_time = 0;
      m_critical_path_wall_clock_time = 0;
      m_critical_path_cpu_time = 0;
      m_num_times_in_critical_path = 0;
      start();
    };
    //Copy constructor
    Timer(const Timer& other) = default;
    Timer& operator=(const Timer& other) = default;
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
      m_last_interval_cpu_time = (static_cast<int64_t>(end_cpu_time.tv_sec - m_begin_cpu_time.tv_sec)*1000000000ll
            + static_cast<int64_t>(end_cpu_time.tv_nsec - m_begin_cpu_time.tv_nsec));
      m_last_interval_wall_clock_time = (static_cast<int64_t>(end_time.tv_sec - m_begin_wall_clock_time.tv_sec)*1000000ll
            + static_cast<int64_t>(end_time.tv_usec - m_begin_wall_clock_time.tv_usec));
      m_cumulative_cpu_time += m_last_interval_cpu_time;
      m_cumulative_wall_clock_time += m_last_interval_wall_clock_time;
    }
    void print_last_interval(const std::string& prefix="", std::ostream& fptr = std::cout) const
    {
      fptr<<"GENOMICSDB_TIMER,";
      if(!prefix.empty())
        fptr << prefix <<",";
      fptr << "Wall-clock time(s),"<< std::setprecision(6) << ((double)m_last_interval_wall_clock_time)/1000000ll << ",Cpu time(s),"
        << ((double)m_last_interval_cpu_time)/1000000000ll << "\n";
    }
    void print(const std::string& prefix="", std::ostream& fptr = std::cout) const
    {
      fptr<<"GENOMICSDB_TIMER,";
      if(!prefix.empty())
        fptr << prefix <<",";
      fptr << "Wall-clock time(s),"<< std::setprecision(6) << ((double)m_cumulative_wall_clock_time)/1000000ll << ",Cpu time(s),"
        << ((double)m_cumulative_cpu_time)/1000000000ll << "\n";
    }
    void print_detail(const std::string& prefix="", std::ostream& fptr = std::cout) const
    {
      fptr<<"GENOMICSDB_TIMER,";
      if(!prefix.empty())
        fptr << prefix <<",";
      fptr << "Wall-clock time(s),"<< std::setprecision(6) << ((double)m_cumulative_wall_clock_time)/1000000ll
        << ",Cpu time(s)," << ((double)m_cumulative_cpu_time)/1000000000ll
        << ",Critical path wall-clock time(s)," << ((double)m_critical_path_wall_clock_time)/1000000ll
        << ",Cpu time(s)," << ((double)m_critical_path_cpu_time)/1000000000ll
        << ",#critical path,"<< m_num_times_in_critical_path << "\n";
    }
    void get_last_interval_times(std::vector<double>& timings, unsigned timer_idx) const
    {
      auto idx = 2u*timer_idx;
      assert(idx + 2u <= timings.size());
      timings[idx] = m_last_interval_cpu_time;
      timings[idx+1] = m_last_interval_wall_clock_time;
    }
    void accumulate(const Timer& other)
    {
      m_cumulative_cpu_time += other.m_cumulative_cpu_time;
      m_cumulative_wall_clock_time += other.m_cumulative_wall_clock_time;
    }
    inline double get_last_interval_wall_clock_time() const { return m_last_interval_wall_clock_time; }
    inline double get_last_interval_cpu_time() const { return m_last_interval_cpu_time; }
    //Critical path updates
    inline void accumulate_critical_path_wall_clock_time(const double val)
    {
      ++m_num_times_in_critical_path;
      m_critical_path_wall_clock_time += val;
    }
    inline void accumulate_critical_path_cpu_time(const double val)
    {
      m_critical_path_cpu_time += val;
    }
  private:
    /*clock_t m_begin_cpu_time;*/
    struct timespec m_begin_cpu_time;
    struct timeval m_begin_wall_clock_time;
    //wall clock in micro-seconds, cpu time in nano-seconds
    uint64_t m_last_interval_cpu_time;
    uint64_t m_last_interval_wall_clock_time;
    uint64_t m_cumulative_cpu_time;
    uint64_t m_cumulative_wall_clock_time;
    //critical path contribution
    uint64_t m_critical_path_wall_clock_time;
    uint64_t m_critical_path_cpu_time;
    uint64_t m_num_times_in_critical_path;
};

struct TimerCompareWallClockTime
{
  bool operator()(const Timer* a, const Timer* b)
  {
    return (a->get_last_interval_wall_clock_time() < b->get_last_interval_wall_clock_time());
  }
  bool operator()(const Timer& a, const Timer& b)
  {
    return (a.get_last_interval_wall_clock_time() < b.get_last_interval_wall_clock_time());
  }
};

#endif
