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

import org.apache.hadoop.io.Writable;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.log4j.Logger;
import org.apache.spark.Partition;
import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.net.InetAddress;

public class GenomicsDBInputSplit extends InputSplit implements Writable {

  // Note that for now the assumption is that there is
  // one GenomicsDB instance per node, hence one host
  // per split
  String[] hosts;
  long length = 0;

  Logger logger = Logger.getLogger(GenomicsDBInputSplit.class);

  /**
   * Default constructor
   */
  public GenomicsDBInputSplit() {
  }

  public GenomicsDBInputSplit(long length) {
    this.length = length;
  }

  public void write(DataOutput dataOutput) throws IOException {
    dataOutput.writeLong(this.length);
  }

  public void readFields(DataInput dataInput) throws IOException {
    hosts = null;
    length = dataInput.readLong();
  }

  public long getLength() throws IOException, InterruptedException {
    return this.length;
  }

  /**
   * Called from {@link org.apache.spark.rdd.NewHadoopRDD#getPreferredLocations(Partition)}
   *
   * @return locations  The values of locations or nodes are passed from host file
   *                   in GenomicsDBConf or hadoopConfiguration in SparkContext
   */
  public String[] getLocations() throws IOException, InterruptedException {
    hosts = new String[1];
    hosts[0] = InetAddress.getLocalHost().getHostName();
    return hosts;
  }
}
