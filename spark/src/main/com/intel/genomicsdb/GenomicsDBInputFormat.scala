/**
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

package com.intel.genomicsdb

import java.util.{List => JavaList}

import genomicsdb.GenomicsDBFeatureReader
import htsjdk.tribble.Feature
import htsjdk.variant.bcf2.BCF2Codec
import org.apache.hadoop.classification.{InterfaceAudience, InterfaceStability}
import org.apache.hadoop.io.LongWritable
import org.apache.hadoop.mapreduce._
import org.apache.spark.TaskContext
import org.apache.spark.annotation.DeveloperApi
import parquet.org.apache.thrift.TException
import parquet.org.apache.thrift.protocol.{TBinaryProtocol, TProtocol}
import parquet.org.apache.thrift.transport.{TSocket, TTransport}

import scala.reflect.ClassTag
import scala.reflect.io.Path

/**
  * For GenomicsDB, the Key.class or K is NULL
  *
  * @tparam VCONTEXT  Variant context type. [[htsjdk.variant.variantcontext.VariantContext]]
  *                   is one of the candidates. See online documentation to know more on this
  */
@InterfaceAudience.Public
@InterfaceStability.Evolving
class GenomicsDBInputFormat[VCONTEXT <: Feature: ClassTag, SOURCE: ClassTag](
    val gConf: GenomicsDBConf)
  extends InputFormat[LongWritable, VCONTEXT] {

  //private val log: Log = LogFactory.getLog(classOf[GenomicsDBInputFormat[VCONTEXT, SOURCE]])
  var recordReader: GenomicsDBRecordReader[VCONTEXT, SOURCE] = _

  /**
    * Lower bound of split size per GenomicsDB partition is 1
    *
    * @return the number of bytes of the minimal split for this format
    */
  protected def getFormatMinSplitSize: Long = 1

  var loaderJsonFile: Path = _
  var queryJsonFile: Path = _
  var hostFile: Path = _

  var featureReader: GenomicsDBFeatureReader[VCONTEXT, SOURCE] = _
  var header: Object = _

  /**
    * Set the loader JSON file if one is provided
    *
    * @param loaderFile  path of the JSON loader file
    */
  def setLoaderJsonFile(loaderFile: Path) = {
    loaderJsonFile = loaderFile
  }

  /**
    * Set the query JSON file if one is provided
    *
    * @param queryFile  path of the JSON query file
    */
  def setQueryJsonFile(queryFile: Path) = {
    queryJsonFile = queryFile
  }

  /**
    * Set the host file if one is provided
    *
    * @param file  full path of the MPI host file
    */
  def setHostFile(file: Path) = {
    hostFile = file
  }

  /**
    * [[getSplits]] method finds the GenomicsDB instances from GenomicsDB configuration
    * and creates one-to-one corresponding input splits. The partitions are in Spark
    * are created based on this number [[com.intel.genomicsdb.GenomicsDBRDD]]
    * Currently, there is 1:1 relation to GenomicsDB data partition per node and
    * input split. This might change in future!!
    *
    * @return  List of input splits
    */
  @DeveloperApi
  @Override
  def getSplits(runningAsService: Boolean = false): List[InputSplit] = {

    // When GenomicsDB runs as a service, get the partition
    // information from the server, else, read from the hosts
    // file provided in the configuration object
    if (runningAsService) {
      // @todo: do nothing yet
      createConnection()
      null
    } else {
      val hostList = gConf.getHosts
      val splits = new Array[InputSplit](hostList.size)

      for (i <- 0 until hostList.size) {
        splits(i) = GenomicsDBInputSplit(Array(hostList.get(i)), gConf)
      }
      splits.toList
    }
  }

  def createRecordReader(
    inputSplit: InputSplit,
    taskContext: TaskContext,
    flag: Int): RecordReader[LongWritable, VCONTEXT] = {

    //val config = taskAttemptContext.getConfiguration

    if (!loaderJsonFile.toString.isEmpty && !queryJsonFile.toString.isEmpty) {
      featureReader = new GenomicsDBFeatureReader(loaderJsonFile.toString,
        queryJsonFile.toString, new BCF2Codec()).asInstanceOf
    } else {
      val tiledbWorkspace = gConf.getTileDBWorkspace
      val arrayName = gConf.getTileDBArrayName
      featureReader = new GenomicsDBFeatureReader(tiledbWorkspace.toString, arrayName,
        gConf.getReferenceGenome, gConf.getVCFHeaderFile, new BCF2Codec()).asInstanceOf
    }
    header = featureReader.getHeader
    recordReader = new GenomicsDBRecordReader[VCONTEXT, SOURCE](inputSplit.asInstanceOf, gConf)
    recordReader
  }

  def getSplits(jobContext: JobContext): JavaList[InputSplit] = {
    // do nothing. never used
    List.empty.asInstanceOf
  }

  /**
    * Return the header from the gVCF from any partition.
    * Header is same across all gVCF files
    *
    * @return
    */
  def getHeader: Object = header


  /**
    * GenomicsDB will be setup as a service to accommodate thousands of compute nodes
    * retreive data from a smaller set of data nodes. In such a scenario, need to create
    * a connection prior to retrieving records and do not assume data will be local
    */
  def createConnection() = {
    try {
      val transport : TTransport = new TSocket(gConf.getHosts.get(0), gConf.getGenomicsDBPort)
      transport.open()

      val protocol: TProtocol = new TBinaryProtocol(transport)
      //val client: GenomicsDBClient = new GenomicsDBService.Client(protocol)
      //perform(client);
      transport.close();
    } catch {
      case e: TException =>
        e.printStackTrace()
    }
  }

  def getRecordReader: GenomicsDBRecordReader[VCONTEXT, SOURCE] = recordReader

  @Override
  def createRecordReader(inputSplit: InputSplit,
                         taskAttemptContext: TaskAttemptContext):
      RecordReader[LongWritable, VCONTEXT] = {
    createRecordReader(inputSplit, null, 1)
  }
}

