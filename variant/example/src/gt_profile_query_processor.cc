#include <iostream>
#include "command_line.h"
#include<math.h>
#include "query_variants.h"
#include "variant_operations.h"

int main(int argc, char** argv) {
  CommandLineOpts cl;
  parse_command_line(argc, argv, cl);
  if(cl.m_workspace == 0 || cl.m_array_name == 0)
  {
    std::cerr << "Missing workspace|array_name\n";
    exit(-1);
  }
  // Create storage manager
  // The input is the path to its workspace (the path must exist).
  VariantStorageManager sm(cl.m_workspace);
  // Create query processor
  VariantQueryProcessor qp(&sm, cl.m_array_name);
  //Setup query config
  VariantQueryConfig query_config;
  query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL"});
  qp.do_query_bookkeeping(qp.get_array_schema(), query_config);
#if 0
  //Iterate over all tiles
  qp.iterate_over_all_tiles(ad_gVCF, query_config);
#endif
  sm.close_array(qp.get_array_descriptor());
  
  return 0;
}
