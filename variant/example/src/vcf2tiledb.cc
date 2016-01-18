#include "vcf2binary.h"
#include "tiledb_loader.h"
#include <mpi.h>

#ifdef USE_GPERFTOOLS
#include "gperftools/profiler.h"
#endif

main(int argc, char** argv)
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
  loader.read_all();
#ifdef USE_GPERFTOOLS
  ProfilerStop();
#endif
  //finalize
  MPI_Finalize();
  return 0;
}
