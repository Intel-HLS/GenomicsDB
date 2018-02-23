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

import com.intel.genomicsdb.reader.GenomicsDBFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.bcf2.BCF2Codec;
import org.apache.hadoop.conf.Configurable;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.mapreduce.*;
import org.apache.log4j.Logger;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class GenomicsDBInputFormat<VCONTEXT extends Feature, SOURCE>
  extends InputFormat<String, VCONTEXT> implements Configurable {

  private GenomicsDBConfiguration genomicsDBConfiguration;
  private Configuration configuration;

  Logger logger = Logger.getLogger(GenomicsDBInputFormat.class);

  /**
   * When this function is called, it is already assumed that configuration
   * object is set
   *
   * @param jobContext  Hadoop Job context passed from newAPIHadoopRDD
   *                    defined in SparkContext
   * @return  Returns a list of input splits
   * @throws FileNotFoundException  Thrown if creaing configuration object fails
   */
  public List<InputSplit> getSplits(JobContext jobContext) throws FileNotFoundException {


    genomicsDBConfiguration = new GenomicsDBConfiguration(configuration);
    genomicsDBConfiguration.setLoaderJsonFile(
      configuration.get(GenomicsDBConfiguration.LOADERJSON));
    genomicsDBConfiguration.setQueryJsonFile(
      configuration.get(GenomicsDBConfiguration.QUERYJSON));
    genomicsDBConfiguration.setHostFile(
      configuration.get(GenomicsDBConfiguration.MPIHOSTFILE));

    List<String> hosts = genomicsDBConfiguration.getHosts();

    ArrayList<InputSplit> inputSplits = new ArrayList<>(hosts.size());
    for (int i = 0; i < hosts.size(); ++i) {
      GenomicsDBInputSplit split = new GenomicsDBInputSplit(hosts.get(i));
      inputSplits.add(split);
    }

    return inputSplits;
  }

  public RecordReader<String, VCONTEXT>
    createRecordReader(InputSplit inputSplit, TaskAttemptContext taskAttemptContext)
      throws IOException, InterruptedException {

    String loaderJson;
    String queryJson;

    GenomicsDBFeatureReader<VCONTEXT, SOURCE> featureReader;
    GenomicsDBRecordReader<VCONTEXT, SOURCE> recordReader;

    if (taskAttemptContext != null) {
      Configuration configuration = taskAttemptContext.getConfiguration();
      loaderJson = configuration.get(GenomicsDBConfiguration.LOADERJSON);
      queryJson = configuration.get(GenomicsDBConfiguration.QUERYJSON);
    } else {
      // If control comes here, means this method is called from
      // GenomicsDBRDD. Hence, the configuration object must be
      // set by setConf method, else this will lead to
      // NullPointerException
      assert(configuration!=null);
      loaderJson = configuration.get(GenomicsDBConfiguration.LOADERJSON);
      queryJson = configuration.get(GenomicsDBConfiguration.QUERYJSON);
    }

    // TODO: fix this once protobuf classes are created
    featureReader = new GenomicsDBFeatureReader<VCONTEXT, SOURCE>(
      loaderJson, "", "", (FeatureCodec<VCONTEXT, SOURCE>) new BCF2Codec());
    recordReader = new GenomicsDBRecordReader<VCONTEXT, SOURCE>(featureReader);
    return recordReader;
  }

  /**
   * default constructor
   */
  public GenomicsDBInputFormat() {
  }

  public GenomicsDBInputFormat(GenomicsDBConfiguration conf) {
    genomicsDBConfiguration = conf;
  }

  /**
   * Set the loader JSON file path
   *
   * @param jsonFile  Full qualified path of the loader JSON file
   * @return  Returns the same object for forward function calls
   */
  public GenomicsDBInputFormat<VCONTEXT, SOURCE> setLoaderJsonFile(String jsonFile) {
    genomicsDBConfiguration.setLoaderJsonFile(jsonFile);
    return this;
  }

  /**
   * Set the query JSON file path
   * @param jsonFile  Full qualified path of the query JSON file
   * @return  Returns the same object for forward function calls
   */
  public GenomicsDBInputFormat<VCONTEXT, SOURCE> setQueryJsonFile(String jsonFile) {
    genomicsDBConfiguration.setQueryJsonFile(jsonFile);
    return this;
  }

  /**
   * Set the host file path
   * @param hostFile  Full qualified path of the hosts file
   * @return  Returns the same object for forward function calls
   * @throws FileNotFoundException thrown if the hosts file is not found
   */
  public GenomicsDBInputFormat<VCONTEXT, SOURCE> setHostFile(String hostFile)
      throws FileNotFoundException {
    genomicsDBConfiguration.setHostFile(hostFile);
    return this;
  }

  @Override
  public void setConf(Configuration configuration) {
    this.configuration = configuration;
  }

  @Override
  public Configuration getConf() {
    return configuration;
  }
}
