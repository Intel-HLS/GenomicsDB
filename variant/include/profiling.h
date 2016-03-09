#ifndef PROFILING_H
#define PROFILING_H

#ifdef DO_PROFILING
#include <vector>
extern uint64_t g_num_disk_loads;
extern uint64_t g_num_cached_loads;
extern uint64_t g_coords_num_disk_loads;
extern uint64_t g_coords_num_cached_loads;
extern uint64_t g_total_num_tiles_loaded;
extern std::vector<uint64_t> g_num_tiles_loaded;
extern std::vector<uint64_t> g_num_segments_loaded;
#endif


#endif
