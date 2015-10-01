#include <iostream>
#include <string>
#include <getopt.h>
#include <mpi.h>
#include "libtiledb_variant.h"

//Enable asserts
#ifdef NDEBUG
#undef NDEBUG
#endif

enum ArgsEnum
{
  ARGS_IDX_SKIP_QUERY_ON_ROOT=1000
};

int main(int argc, char *argv[]) {
  // Define long options
  static struct option long_options[] = 
  {
    {"page-size",1,0,'p'},
    {"output-format",1,0,'O'},
    {"workspace",1,0,'w'},
    {"skip-query-on-root",0,0,ARGS_IDX_SKIP_QUERY_ON_ROOT},
    {"array",1,0,'A'},
    {0,0,0,0},
  };
  int c;
  uint64_t page_size = 0u;
  std::string output_format = "";
  std::string workspace = "";
  std::string array_name = "";
  bool skip_query_on_root = false;
  while((c=getopt_long(argc, argv, "w:A:p:O:", long_options, NULL)) >= 0)
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
      case ARGS_IDX_SKIP_QUERY_ON_ROOT:
        skip_query_on_root = true;
        break;
      default:
        std::cerr << "Unknown command line argument\n";
        exit(-1);
    }
  }
  if(workspace == "" || array_name == "")
  {
    std::cerr << "Missing workspace(-w) or array name (-A)\n";
    return -1;
  }
  if( optind + 2 > argc ) {
    std::cerr << std::endl<< "ERROR: Invalid number of arguments" << std::endl << std::endl;
    std::cout << "Usage: " << argv[0] << " -w <workspace> -A <array name> <start> <end> [ -O <output_format> -p <page_size> ]" << std::endl;
    return -1;
  }
  uint64_t start = std::stoull(std::string(argv[optind]));
  uint64_t end = std::stoull(std::string(argv[optind+1]));
  
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
  /*Create storage manager*/
  StorageManager sm(workspace);
  /*Create query processor*/
  VariantQueryProcessor qp(&sm, array_name);
  //Use VariantQueryConfig to setup query info
  VariantQueryConfig query_config;
  query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "BaseQRankSum", "AD", "PL"});
  query_config.add_column_interval_to_query(start, end);
  qp.do_query_bookkeeping(qp.get_array_schema(), query_config);
  //Variants vector
  std::vector<Variant> variants;
  //Perform query if not root or !skip_query_on_root
  if(my_world_mpi_rank != 0 || !skip_query_on_root)
    qp.gt_get_column_interval(qp.get_array_descriptor(), query_config, 0u, variants);
  //serialized variant data
  std::vector<uint8_t> serialized_buffer;
  serialized_buffer.resize(1000000u);       //1MB, arbitrary value - will be resized if necessary by serialization functions
  uint64_t serialized_length = 0ull;
  for(const auto& variant : variants)
    variant.binary_serialize(serialized_buffer, serialized_length);
  //Gather all serialized lengths (in bytes) at root
  std::vector<uint64_t> lengths_vector(num_mpi_processes, 0ull); 
  assert(MPI_Gather(&serialized_length, 1, MPI_UINT64_T, &(lengths_vector[0]), 1, MPI_UINT64_T, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
  //Total size of gathered data (valid at root only)
  auto total_serialized_size = 0ull;
  //Buffer to receive all gathered data, will be resized at root
  std::vector<uint8_t> receive_buffer(1u);
  //FIXME: MPI uses int for counts and displacements, create multiple batches of gather for large data
  std::vector<int> recvcounts(num_mpi_processes);
  std::vector<int> displs(num_mpi_processes);
  //root
  if(my_world_mpi_rank == 0)
  {
    for(auto val : lengths_vector)
      total_serialized_size += val;
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
  MPI_Gatherv(&(serialized_buffer[0]), serialized_length, MPI_UINT8_T, &(receive_buffer[0]), &(recvcounts[0]), &(displs[0]), MPI_UINT8_T, 0, MPI_COMM_WORLD);
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
    for(const auto& variant : variants)
      variant.print(std::cout, &query_config);
  }
  MPI_Finalize();
  sm.close_array(qp.get_array_descriptor());
  return 0;
}
