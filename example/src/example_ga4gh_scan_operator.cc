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
  // Open arrays in READ mode
  const StorageManager::ArrayDescriptor* ad_gVCF = 
    sm.open_array(cl.m_array_name);

  std::ostream& output_stream = cl.m_output_fstream.is_open() ? cl.m_output_fstream : std::cout;
  // Create query processor
  // The first input is the path to its workspace (the path must exist).
  VariantQueryProcessor qp(cl.m_workspace, sm, ad_gVCF);
#ifdef DO_PROFILING
  sm.m_coords_attribute_idx = qp.get_schema_idx_for_known_field_enum(GVCF_COORDINATES_IDX);
#endif
  //Setup query
  VariantQueryConfig query_config;
  query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL"});
  if(cl.m_position > 0)
    query_config.add_column_interval_to_query(cl.m_position, cl.m_end_position);
  qp.do_query_bookkeeping(ad_gVCF, query_config);
  //Use GA4GH VariantOperators
  GA4GHOperator variant_operator;
  //Do scan and operate
  qp.scan_and_operate(ad_gVCF, query_config, variant_operator, 0u);
  //Print variants - aligned variants (interval splitting, REF, ALT merging etc done)
  for(const auto& variant : variant_operator.get_variants())
    variant.print(output_stream, &query_config);

  sm.close_array(ad_gVCF);
  if(cl.m_output_fstream.is_open())
    cl.m_output_fstream.close();
  if(cl.m_positions_list.is_open())
    cl.m_positions_list.close();

  return 0;
}
