#include <iostream>
#include "command_line.h"
#include "loader.h"
#include "query_processor.h"
#include<math.h>

void GenotypeColumn(QueryProcessor& qp, QueryProcessor::GTProfileStats* stats, const StorageManager::ArrayDescriptor* ad_gVCF,
    uint64_t column, std::ostream& output_stream)
{
  /*Get one column from array*/ 
  QueryProcessor::GTColumn* gt_column = qp.gt_get_column(ad_gVCF, column, stats);
  //Do dummy genotyping operation
  do_dummy_genotyping(gt_column, output_stream);
  delete gt_column;
}

int main(int argc, char** argv) {
  CommandLineOpts cl;
  parse_command_line(argc, argv, cl);
  if(cl.m_workspace == 0 || (cl.m_position == 0ull && !(cl.m_positions_list.is_open()) && !cl.m_do_scan) || cl.m_array_name == 0)
  {
    std::cerr << "Missing workspace|[position or scan]|array_name\n";
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

  std::ostream& output_stream = cl.m_output_fstream.is_open() ? cl.m_output_fstream : std::cout;
  //qp.export_to_CSV(ad_gVCF, "/tmp/blah.csv");
  if(cl.m_do_scan)
    qp.scan_and_operate(ad_gVCF, output_stream);
  else
  {
    QueryProcessor::GTProfileStats stats;
    if(cl.m_positions_list.is_open())
    {
      uint64_t num_queries = 0;
      while(1)
      {
	uint64_t position;
	cl.m_positions_list >> position;
	if(cl.m_positions_list.bad() || cl.m_positions_list.eof() || cl.m_positions_list.fail())
	  break;
	GenotypeColumn(qp, &stats, ad_gVCF, position, output_stream);
	++num_queries;
      }
#ifdef DO_PROFILING
      if(num_queries)
      {
	printf("mean #cells for first sample,mean #cells per sample,mean #cells for exit,mean #tiles per query,,std-dev #cells for first sample,std-dev #cells per sample,std-dev #cells for exit,std-dev #tiles per query\n"); 
	double mean_cells_per_sample = ((double)stats.m_sum_num_cells_touched)/stats.m_num_samples;
	double stddev_cells_per_sample = sqrt(((double)stats.m_sum_sq_num_cells_touched)/stats.m_num_samples - 
	    (mean_cells_per_sample*mean_cells_per_sample));
	double mean_cells_first_sample = ((double)stats.m_sum_num_cells_first_sample)/num_queries;
	double stddev_cells_first_sample = sqrt(((double)stats.m_sum_sq_num_cells_first_sample)/num_queries - 
	    (mean_cells_first_sample*mean_cells_first_sample));
	double mean_cells_last_iter = ((double)stats.m_sum_num_cells_last_iter)/num_queries;
	double stddev_cells_last_iter = sqrt(((double)stats.m_sum_sq_num_cells_last_iter)/num_queries - mean_cells_last_iter*mean_cells_last_iter);
	double mean_tiles_per_query = ((double)stats.m_sum_num_tiles_touched)/num_queries;
	double stddev_tiles_per_query = sqrt(((double)stats.m_sum_sq_num_tiles_touched)/num_queries - (mean_tiles_per_query*mean_tiles_per_query));
	printf("%.2lf,%.2lf,%2.lf,%.2lf,,%.2lf,%.2lf,%2.lf,%.2lf\n",mean_cells_first_sample, mean_cells_per_sample, mean_cells_last_iter, mean_tiles_per_query,stddev_cells_first_sample, stddev_cells_per_sample, stddev_cells_last_iter, stddev_tiles_per_query);
      }
#endif
    }
    else
      GenotypeColumn(qp, &stats, ad_gVCF, cl.m_position, output_stream);
  }

  sm.close_array(ad_gVCF);
  if(cl.m_output_fstream.is_open())
    cl.m_output_fstream.close();
  if(cl.m_positions_list.is_open())
    cl.m_positions_list.close();

  return 0;
}
