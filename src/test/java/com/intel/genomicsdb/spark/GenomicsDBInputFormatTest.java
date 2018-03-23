/*
 * The MIT License (MIT)
 * Copyright 2018 UCLA. License pursuant to original Intel MIT license.
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

import org.testng.Assert;
import org.testng.annotations.Test;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.*;
import org.apache.hadoop.util.ReflectionUtils;

import java.io.File;
import java.net.URL;
import java.net.URISyntaxException;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.InterruptedException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GenomicsDBInputFormatTest {

  @Test(testName = "Testcase0 for creating InputSplits",
      dataProvider = "loaderQueryHostFilesTest0",
      dataProviderClass = GenomicsDBTestUtils.class)
  public void testGetSplits0(String queryPath, String loaderPath, String hostPath) 
              throws IOException, FileNotFoundException, InterruptedException{
    Job job = Job.getInstance();
    Configuration conf = job.getConfiguration();
    conf.set(GenomicsDBConfiguration.LOADERJSON, loaderPath);
    conf.set(GenomicsDBConfiguration.QUERYJSON, queryPath);
    conf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostPath);
    ArrayList<GenomicsDBPartitionInfo> pList = new ArrayList<>(3);
    ArrayList<GenomicsDBQueryInfo> qList = new ArrayList<>(3);
    // test with queries that will not get split up and each
    // query maps to a separate partition
    for(int i=0; i<3; i++) {
      GenomicsDBPartitionInfo p = new GenomicsDBPartitionInfo(i*500000, "hdfs://tmp/ws", "part"+i, "/tmp/test0.vcf.gz");
      GenomicsDBQueryInfo q = new GenomicsDBQueryInfo(i*500000+1000, i*500000+1100);
      pList.add(p);
      qList.add(q);
    }

    GenomicsDBInputFormat format = new GenomicsDBInputFormat();
    format.setConf(conf);
    List<InputSplit> splits = format.getSplits(job);
    Assert.assertEquals(splits.size(), 3);
    for(int i=0; i<splits.size(); i++) {
      GenomicsDBInputSplit gSplit = (GenomicsDBInputSplit)splits.get(i);
      Assert.assertEquals(gSplit.getPartitionInfo(), pList.get(i));
      Assert.assertEquals(gSplit.getQueryInfo(), qList.get(i));
    }
  }

  @Test(testName = "Testcase1 for creating InputSplits",
      dataProvider = "loaderQueryHostFilesTest1",
      dataProviderClass = GenomicsDBTestUtils.class)
  public void testGetSplits1(String queryPath, String loaderPath, String hostPath) 
              throws IOException, FileNotFoundException, InterruptedException{
    Job job = Job.getInstance();
    Configuration conf = job.getConfiguration();
    conf.set(GenomicsDBConfiguration.LOADERJSON, loaderPath);
    conf.set(GenomicsDBConfiguration.QUERYJSON, queryPath);
    conf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostPath);
    ArrayList<GenomicsDBQueryInfo> qList = new ArrayList<>(3);
    // test with queries that will not get split up and all
    // query maps to a single partition
    GenomicsDBPartitionInfo p = new GenomicsDBPartitionInfo(0, "hdfs://tmp/ws", "part", "/tmp/test0.vcf.gz");
    for(int i=0; i<3; i++) {
      GenomicsDBQueryInfo q = new GenomicsDBQueryInfo(i*500000+2000, i*500000+2100);
      qList.add(q);
    }

    GenomicsDBInputFormat format = new GenomicsDBInputFormat();
    format.setConf(conf);
    List<InputSplit> splits = format.getSplits(job);
    Assert.assertEquals(splits.size(), 3);
    for(int i=0; i<splits.size(); i++) {
      GenomicsDBInputSplit gSplit = (GenomicsDBInputSplit)splits.get(i);
      Assert.assertEquals(gSplit.getPartitionInfo(), p);
      Assert.assertEquals(gSplit.getQueryInfo(), qList.get(i));
    }
  }

  @Test(testName = "Testcase2 for creating InputSplits",
      dataProvider = "loaderQueryHostFilesTest2",
      dataProviderClass = GenomicsDBTestUtils.class)
  public void testGetSplits2(String queryPath, String loaderPath, String hostPath) 
              throws IOException, FileNotFoundException, InterruptedException {
    Job job = Job.getInstance();
    Configuration conf = job.getConfiguration();
    conf.set(GenomicsDBConfiguration.LOADERJSON, loaderPath);
    conf.set(GenomicsDBConfiguration.QUERYJSON, queryPath);
    conf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostPath);
    ArrayList<GenomicsDBPartitionInfo> pList = new ArrayList<>(3);
    // test with single query that will get split up and maps to single partition
    GenomicsDBPartitionInfo p = new GenomicsDBPartitionInfo(0, "hdfs://tmp/ws", "part", "/tmp/test0.vcf.gz");
    int qstart = 500;
    int qend = 25000;
    GenomicsDBQueryInfo q = new GenomicsDBQueryInfo(qstart, qend);

    GenomicsDBInputFormat format = new GenomicsDBInputFormat();
    format.setConf(conf);
    List<InputSplit> splits = format.getSplits(job);
    for(int i=0; i<splits.size(); i++) {
      GenomicsDBInputSplit gSplit = (GenomicsDBInputSplit)splits.get(i);
      Assert.assertEquals(gSplit.getPartitionInfo(), p);
      if (i==0) {
        Assert.assertEquals(gSplit.getQueryInfo().getBeginPosition(), qstart);
      }
      else {
        GenomicsDBInputSplit gSplitPrev = (GenomicsDBInputSplit)splits.get(i-1);
        Assert.assertEquals(gSplit.getQueryInfo().getBeginPosition(), 
			gSplitPrev.getQueryInfo().getEndPosition()+1);
      }
      if (i==splits.size()-1) {
        Assert.assertEquals(gSplit.getQueryInfo().getEndPosition(), qend);
      }
    }
  }

  @Test(testName = "Testcase3 for creating InputSplits",
      dataProvider = "loaderQueryHostFilesTest3",
      dataProviderClass = GenomicsDBTestUtils.class)
  public void testGetSplits3(String queryPath, String loaderPath, String hostPath) 
              throws IOException, FileNotFoundException, InterruptedException {
    Job job = Job.getInstance();
    Configuration conf = job.getConfiguration();
    conf.set(GenomicsDBConfiguration.LOADERJSON, loaderPath);
    conf.set(GenomicsDBConfiguration.QUERYJSON, queryPath);
    conf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostPath);
    ArrayList<GenomicsDBPartitionInfo> pList = new ArrayList<>(3);
    // test with single query and non-hdfs compliant data store
    GenomicsDBPartitionInfo p = new GenomicsDBPartitionInfo(0, "/tmp/ws", "part", "/tmp/test0.vcf.gz");
    int qstart = 500;
    int qend = 25000;
    GenomicsDBQueryInfo q = new GenomicsDBQueryInfo(qstart, qend);

    GenomicsDBInputFormat format = new GenomicsDBInputFormat();
    format.setConf(conf);
    List<InputSplit> splits = format.getSplits(job);
    // in this case, number of splits is equal to entries in hostfile
    Assert.assertTrue(splits.size()> 1);
    for(int i=0; i<splits.size(); i++) {
      GenomicsDBInputSplit gSplit = (GenomicsDBInputSplit)splits.get(i);
      // because we're using local data store here, all partitions and
      // queries should be null. GenomicsDB will use the same query json
      // from configuration object for all of them
      Assert.assertEquals(gSplit.getPartitionInfo(), null);
      Assert.assertEquals(gSplit.getQueryInfo(), null);
      // check that locality values have been set according to hostfile
      Assert.assertEquals(gSplit.getLocations()[0], String.valueOf("node"+i));
    }
  }
}
