#include <iostream>
#include "command_line.h"
#include "loader.h"
#include<math.h>
#include "query_variants.h"
#include "variant_operations.h"

int main(int argc, char** argv) {
  CommandLineOpts cl;
  parse_command_line(argc, argv, cl);
  if(cl.m_workspace == 0 || cl.m_array_name == 0)
  {
    std::cerr << "Missing workspace|[position or scan]|array_name\n";
    exit(-1);
  }
  // Create storage manager
  // The input is the path to its workspace (the path must exist).
  StorageManager sm(cl.m_workspace);
  std::ostream& output_stream = cl.m_output_fstream.is_open() ? cl.m_output_fstream : std::cout;
  // Create query processor
  // The first input is the path to its workspace (the path must exist).
  VariantQueryProcessor qp(&sm, cl.m_array_name);
  //Setup query
  VariantQueryConfig query_config;
  query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL", "AD", "GT"});
  if(cl.m_position > 0)
    query_config.add_column_interval_to_query(cl.m_position, cl.m_end_position);
  qp.do_query_bookkeeping(qp.get_array_schema(), query_config);
  //Use GA4GH VariantOperators
  GA4GHOperator variant_operator;
  variant_operator.clear();
#if 0
  //Do scan and operate
  qp.scan_and_operate(ad_gVCF, query_config, variant_operator, 0u);
#endif
  sm.close_array(qp.get_array_descriptor());
  if(cl.m_output_fstream.is_open())
    cl.m_output_fstream.close();
  if(cl.m_positions_list.is_open())
    cl.m_positions_list.close();

  return 0;
}
