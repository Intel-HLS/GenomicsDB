#include <iostream>
#include <string>

#include "libtiledb_variant.h"

int main(int argc, char *argv[]) {
    if( argc < 5 ) {
      char workspace[] = "/mnt/app_hdd/scratch/jagan/TileDB/DB/";
      char array_name[] = "10_DB";
      std::cerr << std::endl<< "ERROR: Invalid number of arguments" << std::endl << std::endl;
      std::cout << "Usage: " << argv[0] << " <workspace> <array name> <start> <end>" << std::endl;
      std::cout << "Example: " << std::endl << "\t" << argv[0]
        << " " << workspace << " " << array_name << " 762588 762600" << std::endl;
      return 0;
    }

    uint64_t start = std::stoull(std::string(argv[3]));
    uint64_t end = std::stoull(std::string(argv[4]));

    bool single_position_queries = false;
    if(argc >= 6 && std::string(argv[5])=="--single-positions")
        single_position_queries = true;
    
    //Use VariantQueryConfig to setup query info
    VariantQueryConfig query_config;
    query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL", "AF", "AN", "AC", "QUAL"});
    if(end > start && single_position_queries)
        //Add interval to query - begin, end
        for( uint64_t i = start; i <= end; ++i )
            query_config.add_column_interval_to_query(i, i);        //single position queries
    else
        query_config.add_column_interval_to_query(start, end);
    
    if(end  == start || single_position_queries)
    {
        Variant variant;
        for( uint64_t i = start; i <= end; ++i ) {
            std::cout << "Position " << i << std::endl;
            db_query_column(argv[1], argv[2], i-start, variant, query_config); 
            std::cout << std::endl;
            variant.print(std::cout, &query_config);
        }
    }
    else
    {
        std::vector<Variant> variants;
        db_query_column_range(argv[1], argv[2], 0ull, variants, query_config);
        for(const auto& variant : variants)
            variant.print(std::cout, &query_config);
    }
    db_cleanup();
}

