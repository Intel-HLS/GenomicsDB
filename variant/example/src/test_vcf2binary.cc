#include "vcf2binary.h"
#include "tiledb_loader.h"

main(int argc, char** argv)
{
  assert(argc >= 2);
  for(auto i=0u;i<1u;++i)
  {
    std::cout << "Column range idx "<<i<<"\n";
    VCF2TileDBLoader loader(argv[1], i);
    loader.read_all();
  }
  return 0;
}
