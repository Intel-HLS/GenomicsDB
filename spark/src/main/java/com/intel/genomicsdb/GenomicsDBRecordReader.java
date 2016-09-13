package com.intel.genomicsdb;

import genomicsdb.GenomicsDBFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.log4j.Logger;

import java.io.IOException;

public class GenomicsDBRecordReader<VCONTEXT extends Feature, SOURCE>
  extends RecordReader<LongWritable, VCONTEXT> {

  private final GenomicsDBFeatureReader<VCONTEXT, SOURCE> featureReader;
  private CloseableTribbleIterator<VCONTEXT> iterator;
  private VCONTEXT currentVariant;
  private int currentKey;

  Logger logger;

  GenomicsDBRecordReader(GenomicsDBFeatureReader<VCONTEXT, SOURCE> featureReader) {
    this.featureReader = featureReader;
    this.currentKey = -1;
    logger = Logger.getLogger(GenomicsDBRecordReader.class);
  }

  public void initialize(InputSplit inputSplit, TaskAttemptContext taskAttemptContext)
    throws IOException, InterruptedException {
    initialize();
  }

  private void initialize() throws IOException {
    logger.info("GenomicsDBRecordReader::initialize() called");

    this.iterator = featureReader.iterator();
    if (iterator.hasNext()) {
      this.currentVariant = iterator.next();
      this.currentKey++;
    }
  }

  public boolean nextKeyValue() throws IOException, InterruptedException {
    if (this.iterator.hasNext()) {
      this.currentVariant = iterator.next();
      this.currentKey++;
      return true;
    } else return false;
  }

  public LongWritable getCurrentKey() throws IOException, InterruptedException {
    return new LongWritable(currentKey);
  }

  public VCONTEXT getCurrentValue() throws IOException, InterruptedException {
    return currentVariant;
  }

  public float getProgress() throws IOException, InterruptedException {
    return 0;
  }

  public void close() throws IOException {
    this.iterator.close();
  }

  public Boolean hasNext() {
    logger.info("GDBRR:hasnext");

    StackTraceElement[] stackTraceElements = Thread.currentThread().getStackTrace();
    for (int i = 0; i < stackTraceElements.length; ++i) {
      logger.info("StackTraceElement[" + i + "]=" + stackTraceElements[i].toString());
    }

    return this.iterator.hasNext();
  }

  public VCONTEXT nextValue() {
    return iterator.next();
  }
}
