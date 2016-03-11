#include <iostream>
#include <string>
#include <getopt.h>

#include "libtiledb_variant.h"

#ifdef USE_GPERFTOOLS
#include "gperftools/profiler.h"
#endif

enum ArgsEnum
{
  ARGS_IDX_TEST_SINGLE_POSITIONS=1000,
  ARGS_IDX_TEST_UPDATE_ROWS,
  ARGS_IDX_TEST_BINARY_SERIALIZATION,
  ARGS_IDX_NUM_ARGS
};

int main(int argc, char *argv[]) {
    // Define long options
    static struct option long_options[] = 
    {
        {"single-positions",0,0,ARGS_IDX_TEST_SINGLE_POSITIONS},
        {"test-update-rows",0,0,ARGS_IDX_TEST_UPDATE_ROWS},
        {"page-size",1,0,'p'},
        {"output-format",1,0,'O'},
        {"test-binary-serialization",0,0,ARGS_IDX_TEST_BINARY_SERIALIZATION},
        {0,0,0,0},
    };
    int c;
    uint64_t page_size = 0u;
    std::string output_format = "";
    bool test_single_positions = false;
    bool test_update_rows = false;
    bool test_binary_serialization = false;
    while((c=getopt_long(argc, argv, "p:O:", long_options, NULL)) >= 0)
    {
        switch(c)
        {
            case 'p':
                page_size = strtoull(optarg, 0, 10);
                break;
            case 'O':
                output_format = std::move(std::string(optarg));
                break;
            case ARGS_IDX_TEST_BINARY_SERIALIZATION:
                test_binary_serialization = true;
                break;
            case ARGS_IDX_TEST_SINGLE_POSITIONS:
                test_single_positions = true;
                break;
            case ARGS_IDX_TEST_UPDATE_ROWS:
                test_update_rows = true;
                break;
            default:
                std::cerr << "Unknown command line argument\n";
                exit(-1);
        }
    }

    if( optind + 4 > argc ) {
      char workspace[] = "/mnt/app_hdd/scratch/jagan/TileDB/DB/";
      char array_name[] = "10_DB";
      std::cerr << std::endl<< "ERROR: Invalid number of arguments" << std::endl << std::endl;
      std::cout << "Usage: " << argv[0] << " <workspace> <array name> <start> <end>" << std::endl;
      std::cout << "Example: " << std::endl << "\t" << argv[0]
        << " " << workspace << " " << array_name << " 762588 762600" << std::endl;
      return 0;
    }
    std::string workspace = argv[optind];
    std::string array_name = argv[optind+1];
    uint64_t start = std::stoull(std::string(argv[optind+2]));
    uint64_t end = std::stoull(std::string(argv[optind+3]));

    //Local TileDB structures
    /*Create storage manager*/
    VariantStorageManager sm(workspace);
    /*Create query processor*/
    VariantQueryProcessor qp(&sm, array_name);

    //Use VariantQueryConfig to setup query info
    VariantQueryConfig query_config;
    //query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL", "GT", "AC", "DP", "PS"});
    query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "BaseQRankSum", "AD", "PL"});
    if(end > start && test_single_positions)
        //Add interval to query - begin, end
        for( uint64_t i = start; i <= end; ++i )
            query_config.add_column_interval_to_query(i, i);        //single position queries
    else
        query_config.add_column_interval_to_query(start, end);
#ifdef USE_GPERFTOOLS
    ProfilerStart("gprofile.log");
#endif
    if(test_single_positions)
    {
        Variant variant;
        for( uint64_t i = start; i <= end; ++i ) {
            std::cout << "Position " << i << std::endl;
            db_query_column(workspace, array_name, i-start, variant, query_config); 
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
          db_query_column_range(workspace, array_name, 0ull, variants, query_config);
          for(const auto& variant : variants)
              variant.print(std::cout, &query_config);
          variants.clear();
          //Bookkeeping must be done BEFORE calling this function
          query_config.update_rows_to_query_to_all_rows();
          std::cout << "Querying all rows\n";
        }
        GA4GHPagingInfo tmp_paging_info;
        GA4GHPagingInfo* paging_info = 0;
        if(page_size)
        {
          paging_info = &tmp_paging_info;
          paging_info->set_page_size(page_size);
        }
        db_query_column_range(workspace, array_name, 0ull, variants, query_config, paging_info);
#ifdef DEBUG
        auto page_idx = 0u;
#endif
        //Multi-page queries
        if(paging_info)
        {
#ifdef DEBUG
          std::cout << "Page "<<page_idx<<" has "<<variants.size()<<" variants\n";
#endif
          std::vector<Variant> tmp_vector;
          while(!(paging_info->is_query_completed()))
          {
            db_query_column_range(workspace, array_name, 0ull, tmp_vector, query_config, paging_info);
            //Move to final vector
            for(auto& v : tmp_vector)
              variants.push_back(std::move(v));
#ifdef DEBUG
            ++page_idx;
            std::cout << "Page "<<page_idx<<" has "<<tmp_vector.size()<<" variants\n";
#endif
            tmp_vector.clear();
          }
        }
        if(test_binary_serialization)
        {
            std::vector<uint8_t> buffer;
            uint64_t offset = 0ull;
            for(const auto& variant : variants)
                variant.binary_serialize(buffer, offset);
            //Clear vector and deserialize into same vector
            auto serialized_length = offset;
            offset = 0ull;
            auto num_variants = variants.size();
            variants.clear();
            variants.resize(num_variants);
            for(auto i=0ull;offset < serialized_length;++i)
                qp.binary_deserialize(variants[i], query_config, buffer, offset);
        }
        print_variants(variants, output_format, query_config, std::cout);
    }
#ifdef USE_GPERFTOOLS
    ProfilerStop();
#endif
    db_cleanup();
    sm.close_array(qp.get_array_descriptor());
    return 0;
}

