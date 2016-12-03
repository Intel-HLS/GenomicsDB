/*
* The MIT License (MIT)
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

package com.intel.genomicsdb;


import com.fasterxml.jackson.databind.ObjectMapper;

import java.io.IOException;
import java.io.Serializable;

/**
 * Maintain global information on how variant data is partitioned
 * across multiple nodes (if loaded with mpiexec)
 *
 * The assumption is that every object described one partition as:
 * {"begin": 0, "workspace":"/home/forall/tiledb-ws",
 *  "array": "test0", "vcf_output_filename":"/home/forall/lihao/mytest/test0.vcf.gz" }
 */
class GenomicsDBPartitionInfo implements Serializable {

  private Integer beginPosition;
  private String tileDBWorkspace;
  private String tileDBArrayName;
  private String VCFOutputFileName;

  public GenomicsDBPartitionInfo(Integer pos,
                                 String workspace, String arrayName, String vcfOutput) {
    beginPosition = pos;
    tileDBWorkspace = workspace;
    tileDBArrayName = arrayName;
    VCFOutputFileName = vcfOutput;
  }

  @Override
  public String toString() {
    String out = beginPosition.toString() + " ";
    out += tileDBWorkspace + " ";
    out += tileDBArrayName + " ";
    out += VCFOutputFileName;
    return out;
  }

  /**
   * Create a JSON string from the partition info object
   *
   * @return  A JSON string of the partition info
   */
  public String toJSON() throws IOException {
    ObjectMapper mapper = new ObjectMapper();
    mapper.writeValueAsString(toString());
    return mapper.toString();
  }
}
