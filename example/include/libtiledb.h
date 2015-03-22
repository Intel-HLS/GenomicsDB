#include "query_processor.h"

class Factory {
  private:
    StorageManager *sm;
    QueryProcessor *qp;
    StorageManager::ArrayDescriptor *ad;

    // path to TileDb workspace
    std::string workspace;
    // Name of the array in the workspace
    std::string array_name;

  public:
    //Stats struct
    QueryProcessor::GTProfileStats stats;

    // Get functions that do lazy initialization of the Tile DB Objects
    StorageManager *getStorageManager(std::string workspace);
    QueryProcessor *getQueryProcessor(std::string workspace);
    StorageManager::ArrayDescriptor *getArrayDescriptor(std::string array_name);
};

extern "C" QueryProcessor::GTColumn *db_query_column(std::string workspace, 
                                                  std::string array_name, 
                                                  uint64_t pos); 
