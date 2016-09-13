package com.intel.genomicsdb;

import genomicsdb.GenomicsDBFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.bcf2.BCF2Codec;
import org.apache.hadoop.conf.Configurable;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.*;
import org.apache.log4j.Logger;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class GenomicsDBInputFormat<VCONTEXT extends Feature, SOURCE>
  extends InputFormat<LongWritable, VCONTEXT> implements Configurable {

  public GenomicsDBConf genomicsDBConf;
  private GenomicsDBFeatureReader<VCONTEXT, SOURCE> featureReader = null;
  private GenomicsDBRecordReader<VCONTEXT, SOURCE> recordReader;

  private Configuration configuration;

  /**
   * When this function is called, it is already assumed that configuration
   * object is set
   *
   * @param jobContext  Hadoop Job context passed from newAPIHadoopRDD
   *                    defined in SparkContext
   * @return
   * @throws IOException
   * @throws InterruptedException
   */
  public List<InputSplit> getSplits(JobContext jobContext)
    throws IOException, InterruptedException {

    System.out.println("Inside GenomicsDBInputFormat::getSplits: started");

    System.out.println("configuration in getsplits: " + configuration.toString());
    System.out.println("configuration in getsplits: " + configuration.get(GenomicsDBConf.LOADERJSON));
    System.out.println("configuration in getsplits: " + configuration.get(GenomicsDBConf.QUERYJSON));
    System.out.println("configuration in getsplits: " + configuration.get(GenomicsDBConf.MPIHOSTFILE));
    genomicsDBConf = new GenomicsDBConf(configuration);
    genomicsDBConf.setLoaderJsonFile(new Path(configuration.get(GenomicsDBConf.LOADERJSON)));
    genomicsDBConf.setQueryJsonFile(new Path(configuration.get(GenomicsDBConf.QUERYJSON)));
    genomicsDBConf.setHostFile(new Path(configuration.get(GenomicsDBConf.MPIHOSTFILE)));
    List<String> hosts = genomicsDBConf.getHosts();

    ArrayList<InputSplit> inputSplits = new ArrayList<InputSplit>(hosts.size());
    for (int i = 0; i < hosts.size(); ++i) {
      GenomicsDBInputSplit split = new GenomicsDBInputSplit();
      split.setLocations(hosts);
      inputSplits.add(split);
    }

    System.out.println("Inside GenomicsDBInputFormat::getSplits: returning splits size: " +
      inputSplits.size());
    return inputSplits;
  }

  public RecordReader<LongWritable, VCONTEXT>
    createRecordReader(InputSplit inputSplit, TaskAttemptContext taskAttemptContext)
      throws IOException, InterruptedException {

    System.out.println("Inside createrecordreader: started");

    if (taskAttemptContext != null) {
      Configuration configuration = taskAttemptContext.getConfiguration();
      String loaderJson = configuration.get(GenomicsDBConf.LOADERJSON);
      String queryJson = configuration.get(GenomicsDBConf.QUERYJSON);

      this.featureReader = new GenomicsDBFeatureReader<VCONTEXT, SOURCE>(
        loaderJson, queryJson, (FeatureCodec<VCONTEXT, SOURCE>) new BCF2Codec());
      this.recordReader = new GenomicsDBRecordReader<VCONTEXT, SOURCE>(this.featureReader);
//      this.recordReader.initialize(inputSplit, taskAttemptContext);
    } else {
      // If control comes here, means this method is called from
      // GenomicsDBRDD. Hence, the configuration object must be
      // set by setConf method, else this will lead to
      // NullPointerException
      assert(configuration!=null);
      String loaderJson = configuration.get(GenomicsDBConf.LOADERJSON);
      String queryJson = configuration.get(GenomicsDBConf.QUERYJSON);

      System.out.println("GenomicsDBInputFormat::createRR::" + loaderJson);
      System.out.println("GenomicsDBInputFormat::createRR::" + queryJson);

      this.featureReader = new GenomicsDBFeatureReader<VCONTEXT, SOURCE>(
        loaderJson, queryJson, (FeatureCodec<VCONTEXT, SOURCE>) new BCF2Codec());
      System.out.println("featurReader dhappa");
      this.recordReader = new GenomicsDBRecordReader<VCONTEXT, SOURCE>(this.featureReader);
    }
    System.out.println("Inside createrecordreader: done");
    return recordReader;
  }

  /**
   * default constructor
   */
  public GenomicsDBInputFormat() {
  }

  public GenomicsDBInputFormat(GenomicsDBConf conf) {
    genomicsDBConf = conf;
  }

  /**
   * Set the loader JSON file path
   *
   * @param jsonFile  Full qualified path of the loader JSON file
   */
  public GenomicsDBInputFormat<VCONTEXT, SOURCE> setLoaderJsonFile(Path jsonFile) {
    genomicsDBConf.setLoaderJsonFile(jsonFile);
    return this;
  }

  public GenomicsDBInputFormat<VCONTEXT, SOURCE> setQueryJsonFile(Path jsonFile) {
    genomicsDBConf.setQueryJsonFile(jsonFile);
    return this;
  }

  public GenomicsDBInputFormat<VCONTEXT, SOURCE> setHostFile(Path hostFile)
      throws FileNotFoundException {
    genomicsDBConf.setHostFile(hostFile);
    return this;
  }

  @Override
  public void setConf(Configuration configuration) {

    Logger logger = Logger.getLogger(GenomicsDBInputFormat.class);
    StackTraceElement[] stackTraceElements = Thread.currentThread().getStackTrace();
    for (int i = 0; i < stackTraceElements.length; ++i) {
      logger.info("StackTraceElement[" + i + "]=" + stackTraceElements[i].toString());
    }

    System.out.println("InputFormat::setConf()");
    String loaderJson = configuration.get(GenomicsDBConf.LOADERJSON);
    String queryJson = configuration.get(GenomicsDBConf.QUERYJSON);
    System.out.println("GenomicsDBInputFormat::setConf::" + loaderJson);
    System.out.println("GenomicsDBInputFormat::setConf::" + queryJson);
    this.configuration = configuration;
  }

  @Override
  public Configuration getConf() {
    return configuration;
  }
}
