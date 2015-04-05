#include "query_variants.h"

class Factory {
  private:
    StorageManager *sm;
    VariantQueryProcessor *qp;
    StorageManager::ArrayDescriptor *ad;

    // path to TileDb workspace
    std::string workspace;
    // Name of the array in the workspace
    std::string array_name;

  public:
    //Stats struct
    GTProfileStats stats;

    // Get functions that do lazy initialization of the Tile DB Objects
    StorageManager *getStorageManager(std::string workspace);
    VariantQueryProcessor *getVariantQueryProcessor(std::string workspace, const StorageManager::ArrayDescriptor* ad);
    StorageManager::ArrayDescriptor *getArrayDescriptor(std::string array_name);
};

extern "C" GTColumn *db_query_column(std::string workspace, 
                                                  std::string array_name, 
                                                  uint64_t pos); 
extern "C" void print_GT_Column(GTColumn *gtc);
