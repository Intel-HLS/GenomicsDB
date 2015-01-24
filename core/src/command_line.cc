#include <iostream>
#include <stdlib.h>
#include <getopt.h>
#include "command_line.h"

void parse_command_line(int argc, char** argv, CommandLineOpts& cl)
{
  static struct option loptions[] =
  {
    {"workspace",0,0,'w'},
    {"csv-file",0,0,'f'},
    {"num-samples",0,0,'N'},
    {"array-name",0,0,'A'},
    {"position",0,0,'p'},
    {0,0,0,0}
  };
  int c = 0;
  while ((c = getopt_long(argc, argv, "w:f:N:p:A:",loptions,NULL)) >= 0) {
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
      case 'p': //position
        cl.m_position = strtoull(optarg, 0, 10);
        break;
      default:
        std::cerr << "Unknown argument "<<argv<<"\n";
        exit(-1);
    }
  }
}
