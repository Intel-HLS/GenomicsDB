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

import org.apache.hadoop.conf.Configuration;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;

/**
 * The configuration class enables users to use Java/Scala
 * to populate the input parameters of GenomicsDB. The first
 * version of the input format are JSON files. Through this
 * class, we provide a programmable interface. Note, that once
 * fully coded, this configuration object can be directly
 * passed from GATK4.0
 */
public class GenomicsDBConfiguration extends Configuration implements Serializable {

  public static final String LOADERJSON = "genomicsdb.input.loaderjsonfile";
  public static final String QUERYJSON = "genomicsdb.input.queryjsonfile";
  public static final String MPIHOSTFILE = "genomicsdb.input.mpi.hostfile";
  public static final String PARTITION_STRATEGY = "genomicsdb.partition.strategy";

  private Boolean produceCombinedVCF = false;
  private Boolean produceTileDBArray = false;
  private Integer segmentSize = 1000;
  private Integer nCellsPerTile = 1000;

  private LinkedList<GenomicsDBPartitionInfo> partitionInfoList = null;

  public GenomicsDBConfiguration(Configuration configuration) throws FileNotFoundException {
    super(configuration);
  }

  /**
   * Constructor with partition information
   * @param configuration  Existing configuration object (can contain Hadoop config values)
   * @param list  Contains partition information
   * @throws FileNotFoundException thrown when loader file not found
   */
  public GenomicsDBConfiguration(Configuration configuration, List<GenomicsDBPartitionInfo> list)
    throws FileNotFoundException {
    super(configuration);
    for (GenomicsDBPartitionInfo info : list) {
      addPartitions(info);
    }
  }

  // <String> left for backward compatibility to Java 7
  private ArrayList<String> hosts = new ArrayList<>();

  public GenomicsDBConfiguration setLoaderJsonFile(String path) {
    set(LOADERJSON, path);
    return this;
  }

  public GenomicsDBConfiguration setQueryJsonFile(String path) {
    set(QUERYJSON, path);
    return this;
  }

  /**
   * Host file contains the hosts where GenomicsDB instances reside.
   * This file can be a replica of the slaves file. For now, we have
   * kept it separate, can be merged later
   *
   * @param path  Full path of the host file
   * @return  GenomicsDBConfiguration object
   * @throws FileNotFoundException  If file not found, throw exception
   */
  public GenomicsDBConfiguration setHostFile(String path) throws FileNotFoundException {
    set(MPIHOSTFILE, path);

    Scanner scanner = new Scanner(new FileInputStream(path));
    while (scanner.hasNextLine()) {
      String host = scanner.nextLine();
      hosts.add(host);
    }
    return this;
  }

  List<String> getHosts() {
    return hosts;
  }

  /**
   * Set the partition strategy for underlying TileDB storage manager
   *
   * @param isRowPartitioned If true, data is partitioned as row major (ROW_MAJOR),
   *                         COL_MAJOR otherwise
   */
  void setPartitionStrategy(Boolean isRowPartitioned) {

    if (isRowPartitioned) {
      set(PARTITION_STRATEGY, String.valueOf(GenomicsDBPartitionStrategy.ROW_MAJOR()));
    } else {
      set(PARTITION_STRATEGY, String.valueOf(GenomicsDBPartitionStrategy.COL_MAJOR()));
    }
  }

  /**
   * Set the partition strategy for underlying TileDB storage manager
   *
   * @param partitionStrategy Partition strategy can be either column major
   *                           (COL_MAJOR) or row major (ROW_MAJOR)
   */
  void setPartitionStrategy(String partitionStrategy) {
    if (partitionStrategy.equals("ROW_MAJOR")) {
      set(PARTITION_STRATEGY, String.valueOf(GenomicsDBPartitionStrategy.ROW_MAJOR()));
    } else if (partitionStrategy.equals("COL_MAJOR")) {
      set(PARTITION_STRATEGY, String.valueOf(GenomicsDBPartitionStrategy.COL_MAJOR()));
    }
  }

  private void addPartitions(GenomicsDBPartitionInfo genomicsDBPartitionInfo) {
    if (partitionInfoList==null) {
      partitionInfoList = new LinkedList<>();
    }
    partitionInfoList.push(genomicsDBPartitionInfo);
  }
}

