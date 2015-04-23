#include <iostream>

#include "libtiledb_variant.h"

Factory f;

StorageManager *Factory::getStorageManager(std::string workspace) {
  if( workspace.compare(this->workspace) != 0 ) {
      if( sm != NULL ) {
          delete sm;
      }
      // Create storage manager
      // The input is the path to its workspace (the path must exist).
      sm = new StorageManager(workspace);
      this->workspace = workspace;
      // Set reset_qp flag since the sm object has changed
      reset_qp = true;
  }
  return sm;
}

VariantQueryProcessor *Factory::getVariantQueryProcessor(std::string workspace, const StorageManager::ArrayDescriptor* ad) {
  if( reset_qp || workspace.compare(this->workspace) != 0 ) {
      // Create query processor
      // The first input is the path to its workspace (the path must exist).
      qp = new VariantQueryProcessor(workspace, *getStorageManager(workspace), ad);
      this->workspace = workspace;
      reset_qp = false;
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

void print_GT_Column(GTColumn *gtc) {
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

extern "C" GTColumn *db_query_column(std::string workspace, 
                                                  std::string array_name, 
                                                  uint64_t pos) {
    // Init Storage Manager object in the Factory class as 
    // both ArrayDescriptor and Query Processor use it 
    StorageManager *sm = f.getStorageManager(workspace);
    StorageManager::ArrayDescriptor *ad = f.getArrayDescriptor(array_name);
    VariantQueryProcessor *qp = f.getVariantQueryProcessor(workspace, ad);

    GTColumn *gtc = qp->gt_get_column(ad, pos, &f.stats);

    return gtc;
}

extern "C" void db_cleanup(std::string workspace, std::string array_name)
{
    StorageManager* sm = f.getStorageManager(workspace);
    StorageManager::ArrayDescriptor *ad = f.getArrayDescriptor(array_name);
    delete f.getVariantQueryProcessor(workspace, ad);
    sm->close_array(f.getArrayDescriptor(array_name));
}
