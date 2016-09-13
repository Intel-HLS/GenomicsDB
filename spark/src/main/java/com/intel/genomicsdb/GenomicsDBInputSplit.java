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
