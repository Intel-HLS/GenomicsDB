#include <iostream>
#include<math.h>
#include "command_line.h"
#include "query_variants.h"
#include "variant_operations.h"

#ifdef USE_GPERFTOOLS
#include "gperftools/profiler.h"
#endif

void GenotypeColumn(VariantQueryProcessor& qp, GTProfileStats* stats, 
    uint64_t column, std::ostream& output_stream, bool do_C_pointer_testing)
{
  VariantQueryConfig query_config;
  query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL"});
  //Add interval to query - begin, end
  query_config.add_column_interval_to_query(column, column);
  qp.do_query_bookkeeping(qp.get_array_schema(), query_config);
  //Variant object
  Variant variant(&query_config);
  variant.resize_based_on_query();
#ifdef DO_PROFILING
  stats->increment_num_queries();
#endif
  /*Get one column from array*/
  qp.gt_get_column(qp.get_array_descriptor(), query_config, 0u, variant, stats);
#ifdef DO_PROFILING
  stats->update_all_from_tmp_counters();
#endif
#if 0
  //Test get_C_pointers() functions
  if(do_C_pointer_testing)
  {
    test_C_pointers(variant);
    std::cout << "Repeating for copy: \n";
    Variant copy;
    copy.copy_from_variant(variant);
    VariantOperations::do_dummy_genotyping(copy, output_stream);
  }
#endif
  //Do dummy genotyping operation
  VariantOperations::do_dummy_genotyping(variant, output_stream);
}

int main(int argc, char** argv) {
  CommandLineOpts cl;
  parse_command_line(argc, argv, cl);
  if(cl.m_workspace == 0 || (cl.m_position == 0ull && !(cl.m_positions_list.is_open()) && !cl.m_do_scan) || cl.m_array_name == 0)
  {
    std::cerr << "Missing workspace|[position or scan]|array_name\n";
    exit(-1);
  }


  std::ostream& output_stream = cl.m_output_fstream.is_open() ? cl.m_output_fstream : std::cout;
  if(cl.m_do_scan)
  {
#if 0
    // Create storage manager
  // The input is the path to its workspace (the path must exist).
  StorageManager sm(cl.m_workspace);
  // Open arrays in READ mode
  const StorageManager::ArrayDescriptor* ad_gVCF = 
    sm.open_array(cl.m_array_name);
    // Create query processor
    // The first input is the path to its workspace (the path must exist).
    VariantQueryProcessor qp(cl.m_workspace, sm, ad_gVCF);
    //Setup query
    VariantQueryConfig query_config;
    query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL"});
    if(cl.m_position > 0)
      query_config.add_column_interval_to_query(cl.m_position, 6000000000ull);  //end is 6B - large value
    qp.do_query_bookkeeping(ad_gVCF, query_config);
    //Use VariantOperators
    DummyGenotypingOperator variant_operator;
    variant_operator.m_output_stream = &(output_stream);
    //Do scan and operate
    qp.scan_and_operate(ad_gVCF, query_config, variant_operator, 0u);
  sm.close_array(ad_gVCF);
#endif
  }
  else
  {
#if 0
    //Compute "optimal" segment size
    const ArraySchema& array_schema = ad_gVCF->array_schema();
    const std::vector<std::pair<double, double> >& dim_domains =
      array_schema.dim_domains();
    uint64_t total_num_samples = dim_domains[0].second - dim_domains[0].first + 1;
    //double expected_num_tiles_per_query =  0.92 + 0.005*total_num_samples;
    //if(expected_num_tiles_per_query < 1)
    //expected_num_tiles_per_query = 1;
#endif
    double expected_num_tiles_per_query = 1;
    uint64_t segment_size = expected_num_tiles_per_query*60000;
    //uint64_t segment_size = 8192;
    if(segment_size > SEGMENT_SIZE)
      segment_size = SEGMENT_SIZE;
    //Create new objects
    StorageManager sm_opt(cl.m_workspace, segment_size);
    // Create query processor
    VariantQueryProcessor qp(&sm_opt, cl.m_array_name);
    //Stats struct
    GTProfileStats stats;
    if(cl.m_positions_list.is_open())
    {
      uint64_t num_queries = 0;
#ifdef USE_GPERFTOOLS
      ProfilerStart("gprofile.log");
#endif
      while(1)
      {
	uint64_t position;
	cl.m_positions_list >> position;
	if(cl.m_positions_list.bad() || cl.m_positions_list.eof() || cl.m_positions_list.fail())
	  break;
	GenotypeColumn(qp, &stats, position, output_stream, cl.m_test_C_pointers);
	++num_queries;
      }
#ifdef USE_GPERFTOOLS
      ProfilerStop();
#endif
    }
    else
      GenotypeColumn(qp, &stats, cl.m_position, output_stream, cl.m_test_C_pointers);
    sm_opt.close_array(qp.get_array_descriptor());
#ifdef DO_PROFILING
    stats.print_stats();
#endif
  }
  if(cl.m_output_fstream.is_open())
    cl.m_output_fstream.close();
  if(cl.m_positions_list.is_open())
    cl.m_positions_list.close();
  return 0;
}
