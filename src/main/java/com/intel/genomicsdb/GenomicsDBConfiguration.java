/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.Iterator;

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

  private ArrayList<GenomicsDBPartitionInfo> partitionInfoList = null;
  private ArrayList<GenomicsDBQueryInfo> queryInfoList = null;
  private long QueryBlockSize = 10000;
  private long QueryBlockSizeMargin = 500;

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
      partitionInfoList = new ArrayList<>();
    }
    partitionInfoList.add(genomicsDBPartitionInfo);
  }

  /**
   * Return partition list; used when creating input splits.
   * @return Returns ArrayList of PartitionInfo objects
   */
  ArrayList<GenomicsDBPartitionInfo> getPartitions() {
    return partitionInfoList;
  }

  /**
   * Return query range list; used when creating input splits.
   * @return Returns ArrayList of QueryRange objects
   */
  ArrayList<GenomicsDBQueryInfo> getQueryRanges() {
    return queryInfoList;
  }

  /**
   * Return value used to determine optimal query size when creating InputSplits
   * @return Returns QueryBlockSize
   */
  long getQueryBlockSize() {
    return QueryBlockSize;
  }

  /**
   * Return value used to determine "slop" for optimal query size when creating InputSplits
   * @return Returns QueryBlockSizeMargin
   */
  long getQueryBlockSizeMargin() {
    return QueryBlockSizeMargin;
  }

  private void readColumnPartitions(JSONObject obj) throws ParseException {
    if (partitionInfoList==null) {
      partitionInfoList = new ArrayList<>();
    }
    JSONArray colPar = (JSONArray)obj.get("column_partitions");
    Iterator<JSONObject> it = colPar.iterator();
    while (it.hasNext()) {
      JSONObject obj0 = (JSONObject)it.next();
      long begin = 0;
      String workspace = null, array = null, vcf_output_filename = null;

      begin = (long)obj0.get("begin");
      workspace = (String)obj0.get("workspace");
      array = (String)obj0.get("array");
      vcf_output_filename = (String)obj0.get("vcf_output_filename");
      partitionInfoList.add(new GenomicsDBPartitionInfo(begin, workspace, array, vcf_output_filename));
    }
  }

  private void readQueryRanges(JSONObject obj) throws ParseException {
    if (queryInfoList==null) {
      queryInfoList = new ArrayList<>();
    }
    // query_column_ranges is a list of lists; we'll only grab first one
    // Assuming here that query file to Spark interface doesn't have a notion
    // of trying to assign certain queries to certain processes or ranks
    assert obj.containsKey("query_column_ranges");
    JSONArray array = (JSONArray)obj.get("query_column_ranges");
    JSONArray firstList = (JSONArray)array.get(0);
    for (Object currElement : firstList) {
      long start = 0, end = 0;
      if (currElement instanceof JSONArray) {
        JSONArray val = (JSONArray)currElement;
	assert val.size() == 2;
        start = (long)val.get(0);
        end = (long)val.get(1);
      }
      else if (currElement instanceof JSONObject) {
        JSONObject val = (JSONObject)currElement;
	assert val.size() == 1;
        start = end = (long)val.get(0);
      }
      else {
        start = end = (long)currElement;
      }
      queryInfoList.add(new GenomicsDBQueryInfo(start, end));
    }

    if(obj.containsKey("query_block_size")) {
      QueryBlockSize = (long)obj.get("query_block_size");
    }
    if(obj.containsKey("query_block_size_margin")) {
      QueryBlockSizeMargin = (long)obj.get("query_block_size_margin");
    }
  }

  /**
   * parse json file to populate ArrayList
   * Assuming here that using the Spark interface implies column partitions
   * @param jsonType json file to use while loading - either LOADERJSON or QUERYJSON
   * @throws FileNotFoundException  Thrown if queryJson file isn't found
   * @throws IOException  Thrown if other IO exception while handling file operations
   * @throws ParseException  Thrown if JSON parsing fails
   */
  void populateListFromJson(String jsonType) 
		  throws FileNotFoundException, IOException, ParseException {
    JSONParser parser = new JSONParser();
    JSONObject obj = (JSONObject)parser.parse(new FileReader(get(jsonType)));

    if (jsonType.equals(LOADERJSON)) {
      readColumnPartitions(obj);
    }
    else if (jsonType.equals(QUERYJSON)) {
      readQueryRanges(obj);
    }
  }
}

