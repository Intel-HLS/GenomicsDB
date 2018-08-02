#include "memory_measure.h"

void read_off_memory_status(statm_t& result, const size_t page_size)
{
  const char* statm_path = "/proc/self/statm";

  FILE *f = fopen(statm_path,"r");
  if(!f){
    perror(statm_path);
    abort();
  }
  if(7 != fscanf(f,"%lu %lu %lu %lu %lu %lu %lu",
        &result.size,&result.resident,&result.share,&result.text,&result.lib,&result.data,&result.dt))
  {
    perror(statm_path);
    abort();
  }
  result.size *= page_size;
  result.resident *= page_size;
  result.share *= page_size;
  result.text *= page_size;
  result.lib *= page_size;
  result.data *= page_size;
  result.dt *= page_size;
  fclose(f);
}
