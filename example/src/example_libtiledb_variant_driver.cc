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
    
    //Use VariantQueryConfig to setup query info
    VariantQueryConfig query_config;
    query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL", "AF", "AN", "AC"});
    //Add interval to query - begin, end
    for( uint64_t i = start; i <= end; ++i )
        query_config.add_column_interval_to_query(i, i);        //currently, only single position queries supported
    
    Variant variant;
    for( uint64_t i = start; i <= end; ++i ) {
        std::cout << "Position " << i << std::endl;
        db_query_column(argv[1], argv[2], i-start, variant, query_config); 
        std::cout << std::endl;
    }
    db_cleanup(argv[1], argv[2]);
}

