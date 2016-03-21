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
#include "command_line.h"
#include<math.h>
#include "query_variants.h"
#include "variant_operations.h"

int main(int argc, char** argv) {
  CommandLineOpts cl;
  parse_command_line(argc, argv, cl);
  if(cl.m_workspace == 0 || cl.m_array_name == 0)
  {
    std::cerr << "Missing workspace|array_name\n";
    exit(-1);
  }
  // Create storage manager
  // The input is the path to its workspace (the path must exist).
  VariantStorageManager sm(cl.m_workspace);
  // Create query processor
  VariantQueryProcessor qp(&sm, cl.m_array_name);
  //Setup query config
  VariantQueryConfig query_config;
  query_config.set_attributes_to_query(std::vector<std::string>{"REF", "ALT", "PL"});
  qp.do_query_bookkeeping(qp.get_array_schema(), query_config);
#if 0
  //Iterate over all tiles
  qp.iterate_over_all_tiles(ad_gVCF, query_config);
#endif
  sm.close_array(qp.get_array_descriptor());
  
  return 0;
}
