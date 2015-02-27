#include <iostream>
#include "command_line.h"
#include "loader.h"
#include "query_processor.h"
#include<math.h>
#include "gt_common.h"

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
  if(cl.m_workspace == 0 || cl.m_array_name == 0)
  {
    std::cerr << "Missing workspace|array_name\n";
    exit(-1);
  }
  // Create storage manager
  // The input is the path to its workspace (the path must exist).
  StorageManager sm(cl.m_workspace);

  // Create query processor
  // The first input is the path to its workspace (the path must exist).
  QueryProcessor qp(cl.m_workspace, sm);

  // Open arrays in READ mode
  const StorageManager::ArrayDescriptor* ad_gVCF = 
    sm.open_array(cl.m_array_name);

  qp.iterate_over_all_cells(ad_gVCF);
#ifdef DO_PROFILING
  printf("COORDS,END,REF,ALT,PL,OFFSETS,NULL\n");
  printf("%.2lf,",((double)g_num_tiles_loaded[GVCF_COORDINATES_IDX])/g_num_segments_loaded[GVCF_COORDINATES_IDX]);
  printf("%.2lf,",((double)g_num_tiles_loaded[GVCF_END_IDX])/g_num_segments_loaded[GVCF_END_IDX]);
  printf("%.2lf,",((double)g_num_tiles_loaded[GVCF_REF_IDX])/g_num_segments_loaded[GVCF_REF_IDX]);
  printf("%.2lf,",((double)g_num_tiles_loaded[GVCF_ALT_IDX])/g_num_segments_loaded[GVCF_ALT_IDX]);
  printf("%.2lf,",((double)g_num_tiles_loaded[GVCF_PL_IDX])/g_num_segments_loaded[GVCF_PL_IDX]);
  printf("%.2lf,",((double)g_num_tiles_loaded[GVCF_OFFSETS_IDX])/g_num_segments_loaded[GVCF_OFFSETS_IDX]);
  printf("%.2lf\n",((double)g_num_tiles_loaded[GVCF_NULL_IDX])/g_num_segments_loaded[GVCF_NULL_IDX]);
#endif
  sm.close_array(ad_gVCF);
  
  return 0;
}
