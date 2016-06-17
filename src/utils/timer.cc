#include "headers.h"
#include "timer.h"

//MacOS does not have a clock_gettime function
//The following function returns wall-clock time and not CPU time
//Any profiling on Mac is pointless for CPU time
#ifdef __MACH__
int clock_gettime(clock_id_t clk_id, struct timespec *tp)
{
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  tp->tv_sec = mts.tv_sec;
  tp->tv_nsec = mts.tv_nsec;
  return 0;
}
#endif
