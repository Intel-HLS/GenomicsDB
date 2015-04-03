#include <iostream>
#include <stdlib.h>
#include <getopt.h>
#include "command_line.h"

enum ArgsIdxEnum
{
  ARGS_IDX_IS_PRESORTED=1000
};

void parse_command_line(int argc, char** argv, CommandLineOpts& cl)
{
  static struct option loptions[] =
  {
    {"array-name",1,0,'A'},
    {"csv-file",1,0,'f'},
    {"num-samples",1,0,'N'},
    {"output",1,0,'o'},
    {"position",1,0,'p'},
    {"position-list",1,0,'P'},
    {"scan",0,0,'S'},
    {"temp-space",1,0,'T'},
    {"workspace",1,0,'w'},
    {"presorted", 0, 0, ARGS_IDX_IS_PRESORTED},
    {0,0,0,0}
  };
  int c = 0;
  char* val;
  while ((c = getopt_long(argc, argv, "A:f:N:o:p:P:Sw:T:",loptions,NULL)) >= 0) {
    switch(c)
    {
      case 'w':
        cl.m_workspace = optarg;
        break;
      case 'f':
        cl.m_csv_filename = optarg;
        break;
      case 'A':
        cl.m_array_name = optarg;
        break;
      case 'N': //#samples
        cl.m_num_samples = strtoull(optarg, 0, 10);
        break;
      case 'o':
        val = optarg;
        cl.m_output_fstream.open(val); 
        break;
      case 'p': //position
        cl.m_position = strtoull(optarg, 0, 10);
        break;
      case 'P':
        val = optarg;
        cl.m_positions_list.open(val); 
        break;
      case 'S':
        cl.m_do_scan = true;
        break;
      case 'T':
        cl.m_temp_space = optarg;
        break;
      case ARGS_IDX_IS_PRESORTED:
	cl.m_is_input_csv_sorted = true;
	break;
      default:
        std::cerr << "Unknown argument "<<argv<<"\n";
        exit(-1);
    }
  }
}
