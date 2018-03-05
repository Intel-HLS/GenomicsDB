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
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import org.apache.hadoop.classification.InterfaceAudience;
import org.apache.hadoop.classification.InterfaceStability;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;

import java.io.IOException;

@InterfaceAudience.Public
@InterfaceStability.Stable
public class GenomicsDBRecordReader<VCONTEXT extends Feature, SOURCE>
  extends RecordReader<String, VCONTEXT>  {

  private final GenomicsDBFeatureReader<VCONTEXT, SOURCE> featureReader;
  private CloseableTribbleIterator<VCONTEXT> iterator;
  private VCONTEXT currentVariant;
  private long currentKey;

  GenomicsDBRecordReader(GenomicsDBFeatureReader<VCONTEXT, SOURCE> featureReader) {
    this.featureReader = featureReader;
    this.currentKey = -1;
  }

  public void initialize(InputSplit inputSplit, TaskAttemptContext taskAttemptContext)
    throws IOException, InterruptedException {
    initialize();
  }

  private void initialize() throws IOException {
    this.iterator = featureReader.iterator();
  }

  public boolean nextKeyValue() throws IOException, InterruptedException {
    if (this.iterator.hasNext()) {
      this.currentVariant = iterator.next();
      this.currentKey++;
      return true;
    } else return false;
  }

  @Override
  public String getCurrentKey() throws IOException, InterruptedException {
    return Long.toString(currentKey);
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
    return this.iterator.hasNext();
  }
}
