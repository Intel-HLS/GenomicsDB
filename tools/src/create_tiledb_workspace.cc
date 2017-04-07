#include <iostream>
#include "tiledb_loader.h"

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Needs 1 argument <workspace_directory>\n";
    exit(-1);
  }
  return VCF2TileDBLoader::create_tiledb_workspace(argv[1]);
}
