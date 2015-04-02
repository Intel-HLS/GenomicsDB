#include <iostream>
#include "libtiledb_variant.h"

int main() {
    for( uint64_t i = 762588; i <= 762600; ++i ) {
        GTColumn *gtc = db_query_column(
                         "/mnt/app_hdd/scratch/jagan/TileDB/DB/", "10_DB", i);
        std::cout << "Position " << i << std::endl;
        print_GT_Column(gtc);
        std::cout << std::endl;
    }
}
