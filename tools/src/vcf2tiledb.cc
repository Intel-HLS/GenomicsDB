/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of 
 * this software and associated documentation files (the "Software"), to deal in 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
 * the Software, and to permit persons to whom the Software is furnished to do so, 
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all 
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "vcf2binary.h"
#include "tiledb_loader.h"
#include <mpi.h>
#include <getopt.h>

#ifdef USE_GPERFTOOLS
#include "gperftools/profiler.h"
#endif

enum VCF2TileDBArgsEnum
{
  VCF2TILEDB_ARG_SPLIT_FILES_IDX=1000,
  VCF2TILEDB_ARG_SPLIT_FILES_PRODUCE_ALL_PARTITIONS_IDX,
  VCF2TILEDB_ARG_SPLIT_FILES_RESULTS_DIRECTORY_IDX,
  VCF2TILEDB_ARG_SPLIT_FILES_SPLIT_OUTPUT_FILENAME_IDX,
  VCF2TILEDB_ARG_SPLIT_FILES_SPLIT_CALLSET_MAPPING_IDX,
  VCF2TILEDB_ARG_VERSION
};

int main(int argc, char** argv)
{
  //Initialize MPI environment
  auto rc = MPI_Init(0, 0);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  //Get my world rank
  int my_world_mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_world_mpi_rank);
  // Define long options
  static struct option long_options[] =
  {
    {"tmp-directory",1,0,'T'},
    {"rank",1,0,'r'},
    {"split-files",0,0,VCF2TILEDB_ARG_SPLIT_FILES_IDX},
    {"split-all-partitions",0,0,VCF2TILEDB_ARG_SPLIT_FILES_PRODUCE_ALL_PARTITIONS_IDX},
    {"split-files-results-directory",1,0,VCF2TILEDB_ARG_SPLIT_FILES_RESULTS_DIRECTORY_IDX},
    {"split-output-filename",1,0,VCF2TILEDB_ARG_SPLIT_FILES_SPLIT_OUTPUT_FILENAME_IDX},
    {"split-callset-mapping-file",0,0,VCF2TILEDB_ARG_SPLIT_FILES_SPLIT_CALLSET_MAPPING_IDX},
    {"version",0,0,VCF2TILEDB_ARG_VERSION},
    {0,0,0,0},
  };
  int c;
  auto split_files = false;
  auto produce_all_partitions = false;
  std::string results_directory;
  std::string split_output_filename;
  auto split_callset_mapping_file = false;
  auto print_version_only = false;
  while((c=getopt_long(argc, argv, "T:r:", long_options, NULL)) >= 0)
  {
    switch(c)
    {
      case 'T':
        g_tmp_scratch_dir = optarg;
        break;
      case 'r':
        my_world_mpi_rank = strtol(optarg, 0, 10);
        break;
      case VCF2TILEDB_ARG_SPLIT_FILES_IDX:
        split_files = true;
        break;
      case VCF2TILEDB_ARG_SPLIT_FILES_PRODUCE_ALL_PARTITIONS_IDX:
        produce_all_partitions = true;
        break;
      case VCF2TILEDB_ARG_SPLIT_FILES_RESULTS_DIRECTORY_IDX:
        results_directory = optarg;
        break;
      case VCF2TILEDB_ARG_SPLIT_FILES_SPLIT_OUTPUT_FILENAME_IDX:
        split_output_filename = optarg;
        break;
      case VCF2TILEDB_ARG_SPLIT_FILES_SPLIT_CALLSET_MAPPING_IDX:
        split_callset_mapping_file = true;
        break;
      case VCF2TILEDB_ARG_VERSION:
        std::cout << GENOMICSDB_VERSION <<"\n";
        print_version_only = true;
        break;
      default:
        std::cerr << "Unknown parameter "<< argv[optind] << "\n";
        exit(-1);
    }
  }
  if(!print_version_only)
  {
    if(optind+1 > argc)
    {
      std::cerr << "Needs <loader_json_config_file>\n";
      exit(-1);
    }
    auto loader_json_config_file = std::move(std::string(argv[optind]));
#ifdef USE_GPERFTOOLS
    ProfilerStart("gprofile.log");
#endif
    //Split files as per the partitions defined - don't load data
    if(split_files)
    {
      GenomicsDBImportConfig loader_config;
      loader_config.read_from_file(loader_json_config_file, my_world_mpi_rank);
      if(loader_config.is_partitioned_by_row())
      {
        std::cerr << "Splitting is available for column partitioning, row partitioning should be trivial if samples are scattered across files. See wiki page https://github.com/Intel-HLS/GenomicsDB/wiki/Dealing-with-multiple-GenomicsDB-partitions for more information\n";
        return 0;
      }
      VidMapper id_mapper = loader_config.get_vid_mapper(); //copy
      //Might specify more VCF files from the command line
      for(auto i=optind+1;i<argc;++i)
        id_mapper.get_or_append_global_file_idx(argv[i]);
      //Single split output
      if(!produce_all_partitions && id_mapper.get_num_files() == 1u && !split_output_filename.empty())
        id_mapper.set_single_split_file_path(0u, split_output_filename);
      std::vector<std::vector<uint8_t>> empty_buffers;
      std::vector<LoaderConverterMessageExchange> empty_exchange;
      const auto& column_partitions = loader_config.get_sorted_column_partitions();
      auto loop_bound = (produce_all_partitions ? column_partitions.size() : 1u);
      for(auto i=0ull;i<loop_bound;++i)
      {
        int rank = produce_all_partitions ? i : my_world_mpi_rank;
        VCF2TileDBConverter converter(loader_config, rank,
            &empty_buffers, &empty_exchange);
        converter.print_all_partitions(results_directory, "", rank);
        if(split_callset_mapping_file)
          id_mapper.write_partition_callsets_json_file(loader_config.get_callset_mapping_file(), results_directory, rank);
      }
      if(split_callset_mapping_file)
        id_mapper.write_partition_loader_json_file(loader_json_config_file, loader_config.get_callset_mapping_file(),
            results_directory, (produce_all_partitions ? column_partitions.size() : 1u), my_world_mpi_rank);
    }
    else
    {
      //Loader object
      VCF2TileDBLoader loader(loader_json_config_file, my_world_mpi_rank);
#ifdef HTSDIR
      loader.read_all();
#endif
    }
#ifdef USE_GPERFTOOLS
    ProfilerStop();
#endif
  }
  //finalize
  MPI_Finalize();
  return 0;
}
