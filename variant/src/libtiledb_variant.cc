#include <iostream>

#include "libtiledb_variant.h"

Factory f;

StorageManager *Factory::getStorageManager(std::string workspace) {
  if( workspace.compare(this->workspace) != 0 ) {
      // Create storage manager
      // The input is the path to its workspace (the path must exist).
      sm = new StorageManager(workspace);
      this->workspace = workspace;
  }
  return sm;
}

QueryProcessor *Factory::getQueryProcessor(std::string workspace) {
  if( workspace.compare(this->workspace) != 0 ) {
      // Create query processor
      // The first input is the path to its workspace (the path must exist).
      qp = new QueryProcessor(workspace, *getStorageManager(workspace));
      this->workspace = workspace;
  }
  return qp;
}

StorageManager::ArrayDescriptor *Factory::getArrayDescriptor(std::string array_name) {
  if( array_name.compare(this->array_name) != 0 ) {
      // Open arrays in READ mode
      ad = sm->open_array(array_name);
      this->array_name = array_name;
  }
  return ad;
}

void print_GT_Column(QueryProcessor::GTColumn *gtc) {
    std::cout << "SAMP \tALT \tREF \tPL \tAF \tAN \tAC" << std::endl;  
    for( int i = 0; i < gtc->ALT_.size(); ++i ) {
        std::cout << i << "\t";
        // Print REF
        std::cout << gtc->REF_[i] << "\t";
        // Print ALT
        for( std::string s : gtc->ALT_[i] ) {
            std::cout << s << ",";
        }
        std::cout << "\t";
        // Print PL
        for( int j : gtc->PL_[i] ) {
            std::cout << j << ",";
        }
        // Print AF
        std::cout << "\t" << gtc->AF_[i];
        // Print AN
        std::cout << "\t" << gtc->AN_[i];
        // Print AC
        std::cout << "\t" << gtc->AC_[i];
        std::cout << std::endl;
    }
}

extern "C" QueryProcessor::GTColumn *db_query_column(std::string workspace, 
                                                  std::string array_name, 
                                                  uint64_t pos) {
    QueryProcessor *qp = f.getQueryProcessor(workspace);
    StorageManager::ArrayDescriptor *ad = f.getArrayDescriptor(array_name);

    QueryProcessor::GTColumn *gtc = qp->gt_get_column(ad, pos, &f.stats);

    return gtc;
}

