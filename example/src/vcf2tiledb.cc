/**
 * The MIT License (MIT)
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

#include "vcf2binary.h"
#include "tiledb_loader.h"
#include <mpi.h>

#ifdef USE_GPERFTOOLS
#include "gperftools/profiler.h"
#endif

int main(int argc, char** argv)
{
  if(argc <= 1)
  {
    std::cerr << "Needs 1 arg <json_config_file>\n";
    exit(-1);
  }
  //Initialize MPI environment
  auto rc = MPI_Init(0, 0);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  //Get number of MPI processes
  int num_mpi_processes = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &num_mpi_processes);
  //Get my world rank
  int my_world_mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_world_mpi_rank);
#ifdef USE_GPERFTOOLS
  ProfilerStart("gprofile.log");
#endif
  //Loader object
  VCF2TileDBLoader loader(argv[1], my_world_mpi_rank);
#ifdef HTSDIR
  loader.read_all();
#endif
#ifdef USE_GPERFTOOLS
  ProfilerStop();
#endif
  //finalize
  MPI_Finalize();
  return 0;
}
