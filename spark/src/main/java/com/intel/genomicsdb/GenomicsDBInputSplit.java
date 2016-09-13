package com.intel.genomicsdb;

import org.apache.hadoop.io.Writable;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.log4j.Logger;
import org.apache.spark.Partition;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.List;


public class GenomicsDBInputSplit extends InputSplit implements Writable {

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
    logger.info("GenomicsDBInputSplit::write() called");
    dataOutput.writeLong(this.length);
  }

  public void readFields(DataInput dataInput) throws IOException {
    // do nothing
    logger.info("GenomicsDBInputSplit::readFields() called");
    hosts = null;
    length = dataInput.readLong();
  }

  public long getLength() throws IOException, InterruptedException {
    return this.length;
  }

  public String[] getLocations() throws IOException, InterruptedException {
    return hosts;
  }

  /**
   * Called from {@link org.apache.spark.rdd.NewHadoopRDD#getPreferredLocations(Partition)}
   *
   * @param locations  The values of locations or nodes are passed from host file
   *                   in GenomicsDBConf or hadoopConfiguration in SparkContext
   */
  public void setLocations(List<String> locations) {
    hosts = new String[locations.size()];
    for (int i = 0; i < locations.size(); ++i)
      hosts[i] = locations.get(i);
  }
}
