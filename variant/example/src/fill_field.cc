#include "query_variants.h"

main(int argc, char** argv)
{
  if(argc < 4)
    assert(0 && "Needs 3 args <workspace> <array_name> <position>");
  StorageManager sm(argv[1]);
  VariantQueryProcessor qp(&sm, argv[2]);
  int64_t position = strtoll(argv[3], 0, 10);
  const ArraySchema* schema_ptr = 0;
  auto status = sm.get_array_schema(qp.get_array_descriptor(), schema_ptr); 
  assert(status == TILEDB_OK);
  //schema_ptr->print();
  VariantQueryConfig query_config;
  query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL"});
  query_config.add_column_interval_to_query(position, position);
  qp.do_query_bookkeeping(*schema_ptr, query_config);
  for(auto i=0u;i<query_config.get_num_queried_attributes();++i)
    std::cout << query_config.get_query_attribute_name(i) << " : "<<query_config.get_schema_idx_for_query_idx(i) << "\n";
  Variant variant(&query_config);
  variant.resize_based_on_query();
  qp.gt_get_column(qp.get_array_descriptor(), query_config, 0, variant);
  variant.print();
  return 0;
}
