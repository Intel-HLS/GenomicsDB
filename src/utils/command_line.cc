/**
 * Copyright (c) 2016 Intel Corporation
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

#include <iostream>
#include <stdlib.h>
#include <getopt.h>
#include "command_line.h"

enum ArgsIdxEnum
{
  ARGS_IDX_IS_PRESORTED=1000,
  ARGS_IDX_TEST_C_POINTERS
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
    {"end-position",1,0,'e'},
    {"position-list",1,0,'P'},
    {"scan",0,0,'S'},
    {"temp-space",1,0,'T'},
    {"workspace",1,0,'w'},
    {"presorted", 0, 0, ARGS_IDX_IS_PRESORTED},
    {"test-c-pointers", 0, 0, ARGS_IDX_TEST_C_POINTERS},
    {0,0,0,0}
  };
  int c = 0;
  char* val;
  while ((c = getopt_long(argc, argv, "A:f:N:o:p:e:P:Sw:T:",loptions,NULL)) >= 0) {
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
      case 'e': //end position
        cl.m_end_position = strtoull(optarg, 0, 10);
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
      case ARGS_IDX_TEST_C_POINTERS:
        cl.m_test_C_pointers = true;
        break;
      default:
        std::cerr << "Unknown argument "<<argv<<"\n";
        exit(-1);
    }
  }
}
