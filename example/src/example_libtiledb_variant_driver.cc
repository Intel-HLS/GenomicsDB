#include <iostream>
#include <string>

#include "libtiledb_variant.h"

int main(int argc, char *argv[]) {
    if( argc < 5 ) {
      char workspace[] = "/mnt/app_hdd/scratch/jagan/TileDB/DB/";
      char array_name[] = "10_DB";
      // char workspace[] = "/mnt/app_hdd/scratch/2arth52g/VCFs/tiledb_csv/arrays/";
      // char array_name[] = "GT10";
      std::cerr << std::endl<< "ERROR: Invalid number of arguments" << std::endl << std::endl;
      std::cout << "Usage: " << argv[0] << " <workspace> <array name> <start> <end>" << std::endl;
      std::cout << "Example: " << std::endl << "\t" << argv[0]
        << " " << workspace << " " << array_name << " 762588 762600" << std::endl;
      return 0;
    }

    uint64_t start = std::stoull(std::string(argv[3]));
    uint64_t end = std::stoull(std::string(argv[4]));
    
    for( uint64_t i = start; i <= end; ++i ) {
        std::cout << "Position " << i << std::endl;
        GTColumn *gtc = db_query_column(argv[1], argv[2], i); 
        print_GT_Column(gtc);
        std::cout << std::endl;
        delete gtc;
    }
    db_cleanup(argv[1], argv[2]);
}

