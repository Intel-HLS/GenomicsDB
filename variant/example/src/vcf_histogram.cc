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
  //Converter object
  VCF2TileDBConverter converter(argv[1], 0);
  converter.create_and_print_histogram(argv[1]);
  return 0;
}
