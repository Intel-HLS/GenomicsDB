#include "query_variants.h"

int main(int argc, char** argv)
{
  if(argc < 4)
    assert(0 && "Needs 3 args <workspace> <array_name> <position> [<end>]");
  StorageManager sm(argv[1]);
  VariantQueryProcessor qp(&sm, argv[2]);
  int64_t position = strtoll(argv[3], 0, 10);
  int64_t end = (argc>=5) ? strtoll(argv[4], 0, 10) : position;
  const ArraySchema* schema_ptr = 0;
  auto status = sm.get_array_schema(qp.get_array_descriptor(), schema_ptr); 
  assert(status == TILEDB_OK);
  //schema_ptr->print();
  VariantQueryConfig query_config;
  query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL"});
  query_config.add_column_interval_to_query(position, end);
  qp.do_query_bookkeeping(*schema_ptr, query_config);
  //Variant variant(&query_config);
  //variant.resize_based_on_query();
  std::vector<Variant> variants;
  variants.clear();
  qp.gt_get_column_interval(qp.get_array_descriptor(), query_config, 0, variants);
  for(auto& var : variants)
  {
    var.set_query_config(&query_config);
    var.print();
  }
  //variant.print();
  return 0;
}
