#include <iostream>
#include <stdlib.h>
#include <getopt.h>
#include "command_line.h"

void parse_command_line(int argc, char** argv, CommandLineOpts& cl)
{
  static struct option loptions[] =
  {
    {"array-name",0,0,'A'},
    {"csv-file",0,0,'f'},
    {"num-samples",0,0,'N'},
    {"output",0,0,'o'},
    {"position",0,0,'p'},
    {"position-list",0,0,'P'},
    {"workspace",0,0,'w'},
    {0,0,0,0}
  };
  int c = 0;
  char* val;
  while ((c = getopt_long(argc, argv, "w:f:N:p:A:o:P:",loptions,NULL)) >= 0) {
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
      default:
        std::cerr << "Unknown argument "<<argv<<"\n";
        exit(-1);
    }
  }
}
