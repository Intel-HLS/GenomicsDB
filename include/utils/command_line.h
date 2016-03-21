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

#ifndef PARSE_COMMAND_LINE_H
#define PARSE_COMMAND_LINE_H

#include<fstream>
class CommandLineOpts
{
  public:
    CommandLineOpts()
    {
      m_do_scan = false;
      m_is_input_csv_sorted = false;
      m_workspace = 0;
      m_csv_filename = 0;
      m_array_name = 0;
      m_num_samples = 0ull;
      m_position = 0ull;
      m_end_position = 600000000000ull; //600B - large number
      m_temp_space = "";
      m_test_C_pointers = false;
    }
    bool m_do_scan;
    bool m_is_input_csv_sorted;
    char* m_workspace;
    char* m_csv_filename;
    char* m_array_name;
    std::ofstream m_output_fstream;
    std::ifstream m_positions_list;
    uint64_t m_num_samples;
    uint64_t m_position;
    uint64_t m_end_position;
    std::string m_temp_space;
    bool m_test_C_pointers;
};

void parse_command_line(int argc, char** argv, CommandLineOpts& cl);

#endif
