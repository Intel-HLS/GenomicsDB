/**
 * @file   example_loader.cc
 * @author Stavros Papadopoulos <stavrosp@csail.mit.edu>
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * Copyright (c) 2014 Stavros Papadopoulos <stavrosp@csail.mit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * This file demonstrates the usage of Loader objects. The loader takes as
 * input a CSV file with raw array data, and creates an array in the
 * native (binary) TileDB format based on the array schema.
 */

#include "loader.h"
#include <iostream>
#include "command_line.h"
int main(int argc, char** argv) {
  CommandLineOpts cl; 
  parse_command_line(argc, argv, cl);
  if(cl.m_workspace == 0 || cl.m_csv_filename == 0 || cl.m_num_samples == 0ull || cl.m_array_name == 0)
  {
    std::cerr << "Missing workspace|csv_filename|num_samples|array_name\n";
    exit(-1);
  }
  try {
    // Create storage manager
    // The input is the path to its workspace (the path must exist).
    StorageManager sm(cl.m_workspace);

    // Create loader
    // The first input is the path to its workspace (the path must exist).
    Loader ld(cl.m_workspace, sm);

    ld.load_CSV_gVCF(cl.m_csv_filename, cl.m_array_name, cl.m_num_samples-1, cl.m_is_input_csv_sorted, cl.m_temp_space);    //argument is max_sample_idx
  } catch(LoaderException& le) {
    std::cout << le.what() << "\n";
  }

  return 0;
}
