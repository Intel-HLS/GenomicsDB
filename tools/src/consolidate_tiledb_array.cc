#include <iostream>
#include "tiledb_loader.h"

int main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Needs 2 arguments <workspace_directory> <array_name>\n";
    exit(-1);
  }
  VCF2TileDBLoader::consolidate_tiledb_array(argv[1], argv[2]);
  return 0;
}
