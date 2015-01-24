#include <iostream>
#include "command_line.h"
#include "loader.h"
#include "query_processor.h"

int main(int argc, char** argv) {
  CommandLineOpts cl;
  parse_command_line(argc, argv, cl);
  if(cl.m_workspace == 0 || cl.m_position == 0ull || cl.m_array_name == 0)
  {
    std::cerr << "Missing workspace|position|array_name\n";
    exit(-1);
  }
  // Create storage manager
  // The input is the path to its workspace (the path must exist).
  StorageManager sm(cl.m_workspace);

  // Create loader
  // The first input is the path to its workspace (the path must exist).
  Loader ld(cl.m_workspace, sm);

  // Create query processor
  // The first input is the path to its workspace (the path must exist).
  QueryProcessor qp(cl.m_workspace, sm);

  // Open arrays in READ mode
  const StorageManager::ArrayDescriptor* ad_gVCF = 
    sm.open_array(cl.m_array_name);

  //qp.export_to_CSV(ad_gVCF, "/tmp/blah.csv");

  //Get one column from array
  QueryProcessor::GTColumn* gt_column = qp.gt_get_column(ad_gVCF, cl.m_position);
  //Do dummy genotyping operation
  do_dummy_genotyping(gt_column);

  sm.close_array(ad_gVCF);
  delete gt_column;

  return 0;
}
