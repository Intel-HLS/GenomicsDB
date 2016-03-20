#include <iostream>
#include <string>
#include <getopt.h>
#include <mpi.h>
#include "libtiledb_variant.h"
#include "json_config.h"
#include "timer.h"
#include "broad_combined_gvcf.h"

#ifdef USE_BIGMPI
#include "bigmpi.h"
#endif

#ifdef USE_GPERFTOOLS
#include "gperftools/profiler.h"
#endif

enum ArgsEnum
{
  ARGS_IDX_SKIP_QUERY_ON_ROOT=1000,
  ARGS_IDX_PRODUCE_BROAD_GVCF,
  ARGS_IDX_PRODUCE_HISTOGRAM,
  ARGS_IDX_DEBUG_PRINT
};

enum CommandsEnum
{
  COMMAND_RANGE_QUERY=0,
  COMMAND_PRODUCE_BROAD_GVCF,
  COMMAND_PRODUCE_HISTOGRAM,
  COMMAND_DEBUG_PRINT
};

#define MegaByte (1024*1024)

#ifdef DO_PROFILING
enum TimerTypesEnum
{
  TIMER_TILEDB_QUERY_RANGE_IDX=0u,
  TIMER_BINARY_SERIALIZATION_IDX,
  TIMER_MPI_GATHER_IDX,
  TIMER_ROOT_BINARY_DESERIALIZATION_IDX,
  TIMER_JSON_PRINTING_IDX,
  TIMER_NUM_TIMERS
};

auto g_timer_names = std::vector<std::string>
{
  "TileDB-query",
  "Binary-serialization",
  "MPI-gather",
  "Binary-deserialization",
  "JSON-printing"
};
#endif

#ifdef NDEBUG
#define ASSERT(X) if(!(X)) { std::cerr << "Assertion failed - exiting\n"; exit(-1); }
#else
#define ASSERT(X) assert(X)
#endif

//id_mapper could be NULL - use for contig/callset name mapping only if non-NULL
void run_range_query(const VariantQueryProcessor& qp, const VariantQueryConfig& query_config, const VidMapper& id_mapper,
    const std::string& output_format, const bool is_partitioned_by_column, int num_mpi_processes, int my_world_mpi_rank, bool skip_query_on_root)
{
  //Check if id_mapper is initialized before using it
  //if(id_mapper.is_initialized())
  GTProfileStats* stats_ptr = 0;
#ifdef DO_PROFILING
  //Performance measurements
  Timer timer;
  std::vector<double> timings(2u*TIMER_NUM_TIMERS, 0);  //[cpu-time,wall-clock time]
  ASSERT(g_timer_names.size() == TIMER_NUM_TIMERS);
  GTProfileStats stats;
  stats_ptr = &stats;
#endif
  //Variants vector
  std::vector<Variant> variants;
  uint64_t num_column_intervals = query_config.get_num_column_intervals();
  std::vector<uint64_t> queried_column_positions(num_column_intervals * 2, 0ull);
  std::vector<uint64_t> query_column_lengths(num_column_intervals, 0ull);
  //Perform query if not root or !skip_query_on_root
  if(my_world_mpi_rank != 0 || !skip_query_on_root)
  {
#ifdef DO_PROFILING
    timer.start();
#endif
    for(auto i=0u;i<query_config.get_num_column_intervals();++i) {
      qp.gt_get_column_interval(qp.get_array_descriptor(), query_config, i, variants, 0, stats_ptr);
      query_column_lengths[i] = variants.size();
      queried_column_positions[i * 2] = query_config.get_column_begin(i);
      queried_column_positions[i * 2 + 1] = query_config.get_column_end(i);
    }

#ifdef DO_PROFILING
    timer.stop();
    timer.get_last_interval_times(timings, TIMER_TILEDB_QUERY_RANGE_IDX);
#endif
  }
#ifdef DO_PROFILING
  timer.start();
#endif
#if VERBOSE>0
  std::cerr << "[Rank "<< my_world_mpi_rank << " ]: Completed query, obtained "<<variants.size()<<" variants\n";
#endif
  //serialized variant data
  std::vector<uint8_t> serialized_buffer;
  serialized_buffer.resize(1000000u);       //1MB, arbitrary value - will be resized if necessary by serialization functions
  uint64_t serialized_length = 0ull;
  for(const auto& variant : variants)
    variant.binary_serialize(serialized_buffer, serialized_length);
#if VERBOSE>0
  std::cerr << "[Rank "<< my_world_mpi_rank << " ]: Completed serialization, serialized data size "
    << std::fixed << std::setprecision(3) << ((double)serialized_length)/MegaByte  << " MBs\n";
#endif
#ifdef DO_PROFILING
  timer.stop();
  timer.get_last_interval_times(timings, TIMER_BINARY_SERIALIZATION_IDX);
  //Gather profiling data at root
  auto num_timing_values_per_mpi_process = 2*(TIMER_BINARY_SERIALIZATION_IDX+1);
  std::vector<double> gathered_timings(num_mpi_processes*num_timing_values_per_mpi_process);
  ASSERT(MPI_Gather(&(timings[0]), num_timing_values_per_mpi_process, MPI_DOUBLE,
        &(gathered_timings[0]), num_timing_values_per_mpi_process, MPI_DOUBLE, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
  timer.start();
#endif
  //Gather all serialized lengths (in bytes) at root
  std::vector<uint64_t> lengths_vector(num_mpi_processes, 0ull); 
  ASSERT(MPI_Gather(&serialized_length, 1, MPI_UNSIGNED_LONG_LONG, &(lengths_vector[0]), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
  //Total size of gathered data (valid at root only)
  auto total_serialized_size = 0ull;
  //Buffer to receive all gathered data, will be resized at root
  std::vector<uint8_t> receive_buffer(1u);
#ifdef USE_BIGMPI
  std::vector<MPI_Count> recvcounts(num_mpi_processes);
  std::vector<MPI_Aint> displs(num_mpi_processes);
#else
  //MPI uses int for counts and displacements
  std::vector<int> recvcounts(num_mpi_processes);
  std::vector<int> displs(num_mpi_processes);
#endif
  //root
  if(my_world_mpi_rank == 0)
  {
    for(auto val : lengths_vector)
      total_serialized_size += val;
#if VERBOSE>0
    std::cerr << "[Rank "<< my_world_mpi_rank << " ]: Gathered lengths, total size "
      << std::fixed << std::setprecision(3) << (((double)total_serialized_size)/MegaByte)<<" MBs\n";
#endif
#ifndef USE_BIGMPI
    if(total_serialized_size >= static_cast<uint64_t>(INT_MAX)) //max 32 bit signed int
    {
      std::cerr << "Serialized size beyond 32-bit int limit - exiting.\n";
      std::cerr <<  "Use the BigMPI library: https://github.com/jeffhammond/BigMPI and recompile the TileDB library with USE_BIGMPI=<path>\n";
      exit(-1);
    }
#endif
    receive_buffer.resize(total_serialized_size);
    auto curr_displ = 0ull;
    //Fill in recvcounts and displs vectors
    for(auto i=0u;i<recvcounts.size();++i)
    {
      auto curr_length = lengths_vector[i];
      recvcounts[i] = curr_length;
      displs[i] = curr_displ;
      curr_displ += curr_length;
    }
  }
  //Gather serialized variant data
#if VERBOSE>0
  if(my_world_mpi_rank == 0)
    std::cerr << "Starting MPIGather at root into buffer of size "<<((double)receive_buffer.size())/MegaByte<<"\n";
#endif
#ifdef USE_BIGMPI
  ASSERT(MPIX_Gatherv_x(&(serialized_buffer[0]), serialized_length, MPI_UNSIGNED_CHAR, &(receive_buffer[0]), &(recvcounts[0]), &(displs[0]), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
#else
  ASSERT(MPI_Gatherv(&(serialized_buffer[0]), serialized_length, MPI_UNSIGNED_CHAR, &(receive_buffer[0]), &(recvcounts[0]), &(displs[0]), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
#endif
  std::vector<uint64_t> gathered_num_column_intervals(num_mpi_processes);
  ASSERT(MPI_Gather(&num_column_intervals, 1, MPI_UNSIGNED_LONG_LONG,
        &(gathered_num_column_intervals[0]), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD) == MPI_SUCCESS);

  uint64_t total_columns = 0ull;
  if(my_world_mpi_rank == 0) {
    for (auto i = 0ull; i < recvcounts.size(); ++i) {
      // Reuse the displs and recvcounts vectors declared for variants serialization
      displs[i] = total_columns;
      total_columns += gathered_num_column_intervals[i];
      recvcounts[i] = gathered_num_column_intervals[i];
    }
  }

  std::vector<uint64_t> gathered_query_column_lengths(total_columns);
#ifdef USE_BIGMPI
  ASSERT(MPIX_Gatherv_x(&(query_column_lengths[0]), num_column_intervals, MPI_UNSIGNED_LONG_LONG,
    &(gathered_query_column_lengths[0]), &(recvcounts[0]), &(displs[0]), MPI_UNSIGNED_LONG_LONG,
    0, MPI_COMM_WORLD) == MPI_SUCCESS);
#else
  ASSERT(MPI_Gatherv(&(query_column_lengths[0]), num_column_intervals, MPI_UNSIGNED_LONG_LONG,
    &(gathered_query_column_lengths[0]), &(recvcounts[0]), &(displs[0]), MPI_UNSIGNED_LONG_LONG,
    0, MPI_COMM_WORLD) == MPI_SUCCESS);
#endif

  // Update the recvcounts and displs but multiply by 2 because
  // each column has the start and end query posistion
  total_columns = 0ull;
  if(my_world_mpi_rank == 0) {
    for (auto i = 0ull; i < recvcounts.size(); ++i) {
      // Reuse the displs and recvcounts vectors declared for variants serialization
      displs[i] = total_columns;
      total_columns += gathered_num_column_intervals[i] * 2;
      recvcounts[i] = gathered_num_column_intervals[i] * 2;
    }
  }

  std::vector<uint64_t> gathered_queried_column_positions(total_columns);
#ifdef USE_BIGMPI
  ASSERT(MPIX_Gatherv_x(&(queried_column_positions[0]), num_column_intervals * 2, MPI_UNSIGNED_LONG_LONG,
    &(gathered_queried_column_positions[0]), &(recvcounts[0]), &(displs[0]), MPI_UNSIGNED_LONG_LONG,
    0, MPI_COMM_WORLD) == MPI_SUCCESS);
#else
  ASSERT(MPI_Gatherv(&(queried_column_positions[0]), num_column_intervals * 2, MPI_UNSIGNED_LONG_LONG,
    &(gathered_queried_column_positions[0]), &(recvcounts[0]), &(displs[0]), MPI_UNSIGNED_LONG_LONG,
    0, MPI_COMM_WORLD) == MPI_SUCCESS);
#endif

#if VERBOSE>0
  if(my_world_mpi_rank == 0)
    std::cerr << "Completed MPI_Gather\n";
#endif
#ifdef DO_PROFILING
  timer.stop();
  timer.get_last_interval_times(timings, TIMER_MPI_GATHER_IDX);
  timer.start();
#endif
  //Deserialize at root
  if(my_world_mpi_rank == 0)
  {
    variants.clear();
    uint64_t offset = 0ull;
    while(offset < total_serialized_size)
    {
      variants.emplace_back();
      auto& variant = variants.back();
      qp.binary_deserialize(variant, query_config, receive_buffer, offset);
    }
#if VERBOSE>0
    std::cerr << "Completed binary deserialization at root\n";
#endif
#ifdef DO_PROFILING
    timer.stop();
    timer.get_last_interval_times(timings, TIMER_ROOT_BINARY_DESERIALIZATION_IDX);
    timer.start();
#endif
    print_variants(variants, output_format, query_config, std::cout, is_partitioned_by_column, &id_mapper,
      gathered_query_column_lengths, gathered_num_column_intervals, gathered_queried_column_positions);
#ifdef DO_PROFILING
    timer.stop();
    timer.get_last_interval_times(timings, TIMER_JSON_PRINTING_IDX);
    std::cerr << "Root received "<< std::fixed << std::setprecision(3) << (((double)total_serialized_size)/MegaByte)
      << " MBs of variant data in binary format\n";
    for(auto i=0u;i<TIMER_NUM_TIMERS;++i)
    {
      std::cerr << g_timer_names[i];
      if(i >= TIMER_MPI_GATHER_IDX) //only root info
        std::cerr << std::fixed << std::setprecision(3) << "," << timings[2*i] << ",," << timings[2*i+1u] << "\n";
      else
      {
        assert(2*i+1 < num_timing_values_per_mpi_process);
        for(auto j=0u;j<gathered_timings.size();j+=num_timing_values_per_mpi_process)
          std::cerr << std::fixed << std::setprecision(3) << "," << gathered_timings[j + 2*i];
        std::cerr << ",";
        for(auto j=0u;j<gathered_timings.size();j+=num_timing_values_per_mpi_process)
          std::cerr << std::fixed << std::setprecision(3) << "," << gathered_timings[j + 2*i + 1];
        std::cerr << "\n";
      }
    }
#endif
  }
}

#if defined(HTSDIR)
void scan_and_produce_Broad_GVCF(const VariantQueryProcessor& qp, const VariantQueryConfig& query_config,
    VCFAdapter& vcf_adapter, const VidMapper& id_mapper,
    int num_mpi_processes, int my_world_mpi_rank, bool skip_query_on_root)
{
  Timer timer;
  BroadCombinedGVCFOperator gvcf_op(vcf_adapter, id_mapper, query_config);
  timer.start();
  qp.scan_and_operate(qp.get_array_descriptor(), query_config, gvcf_op, 0, true);
  timer.stop();
  timer.print("Rank : "+std::to_string(my_world_mpi_rank)+" ", std::cerr);
}
#endif

void debug_print(const VariantQueryProcessor& qp, const VariantQueryConfig& query_config)
{
  VariantCallPrintOperator printer;
  qp.iterate_over_cells(qp.get_array_descriptor(), query_config, printer, 0u);
}

void produce_column_histogram(const VariantQueryProcessor& qp, const VariantQueryConfig& query_config, uint64_t bin_size,
    const std::vector<uint64_t>& num_equi_load_bins)
{
  ColumnHistogramOperator histogram_op(0, 4000000000ull, bin_size);
  qp.iterate_over_cells(qp.get_array_descriptor(), query_config, histogram_op, 0u);
  for(auto val : num_equi_load_bins)
    histogram_op.equi_partition_and_print_bins(val);
}

int main(int argc, char *argv[]) {
  //Initialize MPI environment
  auto rc = MPI_Init(0, 0);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  //Get number of MPI processes
  int num_mpi_processes = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &num_mpi_processes);
  //Get my world rank
  int my_world_mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_world_mpi_rank);
#ifdef DEBUG
  //Print host, rank and LD_LIBRARY_PATH
  std::vector<char> hostname;
  hostname.resize(500u);
  gethostname(&(hostname[0]), hostname.size());
  if(my_world_mpi_rank == 0)
    std::cerr << "#processes "<<num_mpi_processes<<"\n";
  std::cerr << "Host : "<< &(hostname[0]) << " rank "<< my_world_mpi_rank << " LD_LIBRARY_PATH= "<<getenv("LD_LIBRARY_PATH") << "\n";
#endif
  // Define long options
  static struct option long_options[] = 
  {
    {"page-size",1,0,'p'},
    {"output-format",1,0,'O'},
    {"workspace",1,0,'w'},
    {"json-config",1,0,'j'},
    {"loader-json-config",1,0,'l'},
    {"segment-size",1,0,'s'},
    {"skip-query-on-root",0,0,ARGS_IDX_SKIP_QUERY_ON_ROOT},
    {"produce-Broad-GVCF",0,0,ARGS_IDX_PRODUCE_BROAD_GVCF},
    {"produce-histogram",0,0,ARGS_IDX_PRODUCE_HISTOGRAM},
    {"debug-print",0,0,ARGS_IDX_DEBUG_PRINT},
    {"array",1,0,'A'},
    {0,0,0,0},
  };
  int c;
  uint64_t page_size = 0u;
  std::string output_format = "";
  std::string workspace = "";
  std::string array_name = "";
  std::string json_config_file = "";
  std::string loader_json_config_file = "";
  bool skip_query_on_root = false;
  unsigned command_idx = COMMAND_RANGE_QUERY;
  size_t segment_size = 10u*1024u*1024u; //in bytes = 10MB
  while((c=getopt_long(argc, argv, "j:l:w:A:p:O:s:", long_options, NULL)) >= 0)
  {
    switch(c)
    {
      case 'p':
        page_size = strtoull(optarg, 0, 10);
        std::cerr << "WARNING: page size is ignored for now\n";
        break;
      case 'O':
        output_format = std::move(std::string(optarg));
        break;
      case 'w':
        workspace = std::move(std::string(optarg));
        break;
      case 'A':
        array_name = std::move(std::string(optarg));
        break;
      case 's':
        segment_size = strtoull(optarg, 0, 10);
        break;
      case ARGS_IDX_SKIP_QUERY_ON_ROOT:
        skip_query_on_root = true;
        break;
      case ARGS_IDX_PRODUCE_BROAD_GVCF:
        command_idx = COMMAND_PRODUCE_BROAD_GVCF;
        break;
      case ARGS_IDX_PRODUCE_HISTOGRAM:
        command_idx = COMMAND_PRODUCE_HISTOGRAM;
        break;
      case 'j':
        json_config_file = std::move(std::string(optarg));
        break;
      case ARGS_IDX_DEBUG_PRINT:
        command_idx = COMMAND_DEBUG_PRINT;
        break;
      case 'l':
        loader_json_config_file = std::move(std::string(optarg));
        break;
      default:
        std::cerr << "Unknown command line argument\n";
        exit(-1);
    }
  }
  //Use VariantQueryConfig to setup query info
  VariantQueryConfig query_config;
  //Vid mapping
  FileBasedVidMapper id_mapper;
  JSONLoaderConfig loader_config;
  auto file_id_mapper_ptr = &id_mapper;
#ifdef HTSDIR
  VCFAdapter vcf_adapter;
#endif
  //If JSON file specified, read workspace, array_name, rows/columns/fields to query from JSON file
  if(json_config_file != "")
  {
    JSONBasicQueryConfig* json_config_ptr = 0;
    JSONBasicQueryConfig range_query_config;
#if defined(HTSDIR)
    JSONVCFAdapterQueryConfig scan_config;
#endif
    //If loader JSON passed as argument, initialize id_mapper from this file
    if(loader_json_config_file.length())
    {
      loader_config.read_from_file(loader_json_config_file, &id_mapper, my_world_mpi_rank);
      //Do not initialize FileBasedVidMapper from the query JSON
      file_id_mapper_ptr = 0;
    }
    switch(command_idx)
    {
      case COMMAND_PRODUCE_BROAD_GVCF:
#if defined(HTSDIR)
        scan_config.read_from_file(json_config_file, query_config, vcf_adapter, file_id_mapper_ptr, output_format, my_world_mpi_rank);
        json_config_ptr = static_cast<JSONBasicQueryConfig*>(&scan_config);
#else
        std::cerr << "Cannot produce Broad's combined GVCF without htslib. Re-compile with HTSDIR variable set\n";
        exit(-1);
#endif
        break;
      default:
        range_query_config.read_from_file(json_config_file, query_config, &id_mapper, my_world_mpi_rank, &loader_config);
        json_config_ptr = &range_query_config;
        break;
    }
    ASSERT(json_config_ptr);
    workspace = json_config_ptr->get_workspace(my_world_mpi_rank);
    array_name = json_config_ptr->get_array_name(my_world_mpi_rank);
  }
  else
  {
    if( optind + 2 > argc ) {
      std::cerr << std::endl<< "ERROR: Invalid number of arguments" << std::endl << std::endl;
      std::cout << "Usage: " << argv[0] << "  ( -j <json_config_file> | -w <workspace> -A <array name> <start> <end> ) [ -O <output_format> -p <page_size> ]" << std::endl;
      return -1;
    }
    uint64_t start = std::stoull(std::string(argv[optind]));
    uint64_t end = std::stoull(std::string(argv[optind+1])); 
    query_config.add_column_interval_to_query(start, end);
    switch(command_idx)
    {
      case COMMAND_PRODUCE_BROAD_GVCF:
        std::cerr << "To produce Broad's combined GVCF, you need to pass parameters through a JSON file, exiting\n";
        exit(-1);
        break;
      case COMMAND_PRODUCE_HISTOGRAM:
        break;  //no attributes
      case COMMAND_DEBUG_PRINT:
        query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT"});
        break;
      default:
        query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "BaseQRankSum", "AD", "PL"});
        break;
    }
  }
  if(workspace == "" || array_name == "")
  {
    std::cerr << "Missing workspace(-w) or array name (-A)\n";
    return -1;
  }
#ifdef USE_GPERFTOOLS
  ProfilerStart("gprofile.log");
#endif
#if VERBOSE>0
  std::cerr << "Segment size: "<<segment_size<<" bytes\n";
#endif
  /*Create storage manager*/
  VariantStorageManager sm(workspace, segment_size);
  /*Create query processor*/
  VariantQueryProcessor qp(&sm, array_name);
  qp.do_query_bookkeeping(qp.get_array_schema(), query_config);
  switch(command_idx)
  {
    case COMMAND_RANGE_QUERY:
      run_range_query(qp, query_config, static_cast<const VidMapper&>(id_mapper), output_format,
          loader_config.is_partitioned_by_column(), num_mpi_processes, my_world_mpi_rank, skip_query_on_root);
      break;
    case COMMAND_PRODUCE_BROAD_GVCF:
#if defined(HTSDIR)
      scan_and_produce_Broad_GVCF(qp, query_config, vcf_adapter, static_cast<const VidMapper&>(id_mapper),
          num_mpi_processes, my_world_mpi_rank, skip_query_on_root);
#endif
      break;
    case COMMAND_PRODUCE_HISTOGRAM:
      produce_column_histogram(qp, query_config, 100, std::vector<uint64_t>({ 128, 64, 32, 16, 8, 4, 2 }));
      break;
    case COMMAND_DEBUG_PRINT:
      debug_print(qp, query_config);
      break;
  }
#ifdef USE_GPERFTOOLS
  ProfilerStop();
#endif

  MPI_Finalize();
  sm.close_array(qp.get_array_descriptor());
  return 0;
}
