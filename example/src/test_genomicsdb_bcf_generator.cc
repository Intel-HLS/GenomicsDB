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

#include <getopt.h>
#include "headers.h"
#include "genomicsdb_bcf_generator.h"
#include <mpi.h>

int main(int argc, char *argv[]) {
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
    {"page-size",1,0,'p'},
    {"rank",1,0,'r'},
    {"output-format",1,0,'O'},
    {"json-config",1,0,'j'},
    {"loader-json-config",1,0,'l'},
    {0,0,0,0},
  };
  int c;
  uint64_t page_size = 0u;
  std::string output_format = "";
  std::string json_config_file = "";
  std::string loader_json_config_file = "";
  while((c=getopt_long(argc, argv, "j:l:p:r:O:", long_options, NULL)) >= 0)
  {
    switch(c)
    {
      case 'p':
        page_size = strtoull(optarg, 0, 10);
        break;
      case 'r':
        my_world_mpi_rank = strtoull(optarg, 0 ,10);
        break;
      case 'O':
        output_format = std::move(std::string(optarg));
        break;
      case 'j':
        json_config_file = std::move(std::string(optarg));
        break;
      case 'l':
        loader_json_config_file = std::move(std::string(optarg));
        break;
      default:
        std::cerr << "Unknown command line argument\n";
        exit(-1);
    }
  }
  std::vector<uint8_t> buffer(page_size > 0u ? page_size : 100u); 
  //assert(json_config_file.length() > 0u && loader_json_config_file.length() > 0u);
  GenomicsDBBCFGenerator bcf_reader(loader_json_config_file, json_config_file, my_world_mpi_rank, page_size, std::max<size_t>(page_size, 1024u),
      output_format.c_str());
  while(!(bcf_reader.end()))
  {
    auto num_bytes_read = bcf_reader.read_and_advance(&(buffer[0]), 0u, buffer.size());
    if(num_bytes_read > 0u)
      fwrite(&(buffer[0]), 1u, num_bytes_read, stdout);
  }
  MPI_Finalize();
  return 0;
}
