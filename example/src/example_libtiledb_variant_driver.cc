#include <iostream>
#include "libtiledb_variant.h"

int main() {
    for( uint64_t i = 762588; i <= 762600; ++i ) {
        std::cout << "Position " << i << std::endl;
        QueryProcessor::GTColumn *gtc = db_query_column(
                         "/mnt/app_hdd/scratch/jagan/TileDB/DB/", "2_ICGC", i);
        print_GT_Column(gtc);
        std::cout << std::endl;
    }
    for( uint64_t i = 762344; i <= 762347; ++i ) {
        std::cout << "Position " << i << std::endl;
        QueryProcessor::GTColumn *gtc = db_query_column(
                         "/mnt/app_hdd/scratch/jagan/TileDB/DB/", "2_ICGC", i);
        print_GT_Column(gtc);
        std::cout << std::endl;
    }
}
