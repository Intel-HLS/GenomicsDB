#include <iostream>
#include <stdlib.h>
#include <getopt.h>

void parse_command_line(int argc, char** argv, char*& workspace, char*& csv_filename, uint64_t& value)
{
  static struct option loptions[] =
  {
    {"workspace",0,0,'w'},
    {"csv-file",0,0,'f'},
    {"num-samples",0,0,'N'},
    {"position",0,0,'p'},
    {0,0,0,0}
  };
  int c = 0;
  while ((c = getopt_long(argc, argv, "w:f:N:p:",loptions,NULL)) >= 0) {
    switch(c)
    {
      case 'w':
        workspace = optarg;
        break;
      case 'f':
        csv_filename = optarg;
        break;
      case 'N': //#samples
      case 'p': //position
        value = strtoull(optarg, 0, 10);
        break;
      default:
        std::cerr << "Unknown argument "<<argv<<"\n";
        exit(-1);
    }
  }
}
