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

extern "C" void db_query_column(std::string workspace, std::string array_name, 
        uint64_t query_interval_idx, Variant& variant, VariantQueryConfig& query_config) {
    // Init Storage Manager object in the Factory class as 
    // both ArrayDescriptor and Query Processor use it 
    StorageManager *sm = f.getStorageManager(workspace);
    StorageManager::ArrayDescriptor *ad = f.getArrayDescriptor(array_name);
    VariantQueryProcessor *qp = f.getVariantQueryProcessor(workspace, ad);
    //Do book-keeping, if not already done
    if(!query_config.is_bookkeeping_done())
    {
      qp->do_query_bookkeeping(ad, query_config);
      variant = std::move(Variant(&query_config));
      variant.resize_based_on_query();
    }
    qp->gt_get_column(ad, query_config, query_interval_idx, variant, &f.stats);
    variant.print(std::cout);
}

extern "C" void db_query_column_range(std::string workspace, std::string array_name, 
        uint64_t query_interval_idx, std::vector<Variant>& variants, VariantQueryConfig& query_config) {
    // Init Storage Manager object in the Factory class as 
    // both ArrayDescriptor and Query Processor use it 
    StorageManager *sm = f.getStorageManager(workspace);
    StorageManager::ArrayDescriptor *ad = f.getArrayDescriptor(array_name);
    VariantQueryProcessor *qp = f.getVariantQueryProcessor(workspace, ad);
    //Do book-keeping, if not already done
    if(!query_config.is_bookkeeping_done())
        qp->do_query_bookkeeping(ad, query_config);
    qp->gt_get_column_interval(ad, query_config, query_interval_idx, variants, &f.stats);
    for(const auto& variant : variants)
        variant.print(std::cout, &query_config);
}
extern "C" void db_cleanup(std::string workspace, std::string array_name)
{
    StorageManager* sm = f.getStorageManager(workspace);
    StorageManager::ArrayDescriptor *ad = f.getArrayDescriptor(array_name);
    delete f.getVariantQueryProcessor(workspace, ad);
    sm->close_array(f.getArrayDescriptor(array_name));
}
