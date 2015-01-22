#include <iostream>
#include "command_line.h"
#include "loader.h"
#include "query_processor.h"

int main(int argc, char** argv) {
  char* workspace = 0;
  char* csv_filename = 0;
  uint64_t position = 0ull;
  parse_command_line(argc, argv, workspace, csv_filename, position);
  if(workspace == 0 || csv_filename == 0 || position == 0ull)
  {
    std::cerr << "Missing workspace|csv_filename|position\n";
    exit(-1);
  }
  // Create storage manager
  // The input is the path to its workspace (the path must exist).
  StorageManager sm(workspace);

  // Create loader
  // The first input is the path to its workspace (the path must exist).
  Loader ld(workspace, sm);

  // Create query processor
  // The first input is the path to its workspace (the path must exist).
  QueryProcessor qp(workspace, sm);

  // Open arrays in READ mode
  const StorageManager::ArrayDescriptor* ad_gVCF = 
    sm.open_array("GEN");     //hard coded to the value in loader

  //Get one column from array
  QueryProcessor::GTColumn* gt_column = qp.gt_get_column(ad_gVCF, position);
  //Do dummy genotyping operation
  do_dummy_genotyping(gt_column);

  sm.close_array(ad_gVCF);

  return 0;
}
