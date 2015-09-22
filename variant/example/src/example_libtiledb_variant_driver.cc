#include <iostream>
#include <string>

#include "libtiledb_variant.h"

#ifdef USE_GPERFTOOLS
#include "gperftools/profiler.h"
#endif

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
    bool test_update_rows = false;
    if(argc >= 6) 
      if(std::string(argv[5])=="--single-positions")
        single_position_queries = true;
      else
        if(std::string(argv[5])=="--test-update-rows")
          test_update_rows = true;
    
    //Use VariantQueryConfig to setup query info
    VariantQueryConfig query_config;
    //query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL", "GT", "AC", "DP", "PS"});
    query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL"});
    if(end > start && single_position_queries)
        //Add interval to query - begin, end
        for( uint64_t i = start; i <= end; ++i )
            query_config.add_column_interval_to_query(i, i);        //single position queries
    else
        query_config.add_column_interval_to_query(start, end);
#ifdef USE_GPERFTOOLS
    ProfilerStart("gprofile.log");
#endif
    if(single_position_queries)
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
        if(test_update_rows)
        {
          std::cout << "Querying only row 0\n";
          query_config.set_rows_to_query(std::vector<int64_t>(1u, 0ll));     //query row 0 only
          db_query_column_range(argv[1], argv[2], 0ull, variants, query_config);
          for(const auto& variant : variants)
              variant.print(std::cout, &query_config);
          variants.clear();
          //Bookkeeping must be done BEFORE calling this function
          query_config.update_rows_to_query_to_all_rows();
          std::cout << "Querying all rows\n";
        }
        db_query_column_range(argv[1], argv[2], 0ull, variants, query_config);
        for(const auto& variant : variants)
            variant.print(std::cout, &query_config);
    }
#ifdef USE_GPERFTOOLS
    ProfilerStop();
#endif
    db_cleanup();
}

