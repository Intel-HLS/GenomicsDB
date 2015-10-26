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
      m_begin_cpu_time = clock();
      gettimeofday(&m_begin_wall_clock_time, 0);
    }
    inline void stop()
    {
      m_last_interval_cpu_time = ((double)(clock() - m_begin_cpu_time))/CLOCKS_PER_SEC;
      struct timeval end_time;
      gettimeofday(&end_time, 0);
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
    void get_last_interval_times(std::vector<double>& timings, unsigned timer_idx) const
    {
      auto idx = 2u*timer_idx;
      assert(idx + 2u <= timings.size());
      timings[idx] = m_last_interval_cpu_time;
      timings[idx+1] = m_last_interval_wall_clock_time;
    }
  private:
    clock_t m_begin_cpu_time;
    struct timeval m_begin_wall_clock_time;
    //seconds
    double m_last_interval_cpu_time;
    double m_last_interval_wall_clock_time;
    double m_cumulative_cpu_time;
    double m_cumulative_wall_clock_time;
};

#endif
