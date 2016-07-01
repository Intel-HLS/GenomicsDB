#ifndef MEMORY_MEASURE_H
#define MEMORY_MEASURE_H

#include <stdio.h>
#include <stdlib.h>

typedef struct {
  unsigned long size,resident,share,text,lib,data,dt;
} statm_t;

void read_off_memory_status(statm_t& result, const size_t page_size=4096u);

#endif
