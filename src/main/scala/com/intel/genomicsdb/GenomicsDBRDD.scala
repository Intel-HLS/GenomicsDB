/**
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

package com.intel.genomicsdb

import com.intel.genomicsdb.spark.{GenomicsDBConfiguration, GenomicsDBInputFormat, GenomicsDBRecordReader}
import com.intel.genomicsdb.spark.GenomicsDBInputFormat
import htsjdk.tribble.Feature
import org.apache.hadoop.conf.Configuration
import org.apache.hadoop.io.Writable
import org.apache.hadoop.mapreduce.{InputSplit, RecordReader}
import org.apache.spark._
import org.apache.spark.annotation.DeveloperApi
import org.apache.spark.api.java.JavaRDD
import org.apache.spark.broadcast.Broadcast
import org.apache.spark.rdd.RDD
import org.apache.spark.storage.StorageLevel

import scala.annotation.meta.param
import scala.reflect.ClassTag

/**
  * Resilient distributed structures containing variant data. The
  * variants are populated from GenomicsDB (data source)
  *
  * @param gc GenomicsDB context contains the Spark Context
  * @param parallelism  Number of RDD partitions
  */
class GenomicsDBRDD[VCONTEXT <: Feature: ClassTag, SOURCE: ClassTag](
    @(transient @param) gc: GenomicsDBContext,
    parallelism: Int,
    runningAsService: Boolean = false)
  extends RDD[VCONTEXT](gc.sparkContext, Nil) {

  private val shouldCloneJobConf = sparkContext.getConf.getBoolean("spark.hadoop.cloneConf",
    defaultValue = false)

  private val confBroadcast = gc.sparkContext.broadcast(
    new GenomicsDBConfiguration(gc.sparkContext.hadoopConfiguration))

  val loaderJsonBroadcast: Broadcast[String] = sparkContext.broadcast(
    confBroadcast.value.get(GenomicsDBConfiguration.LOADERJSON))
  val queryJsonBroadcast: Broadcast[String] = sparkContext.broadcast(
    confBroadcast.value.get(GenomicsDBConfiguration.QUERYJSON))
  val hostsBroadcast: Broadcast[String] = sparkContext.broadcast(
    confBroadcast.value.get(GenomicsDBConfiguration.MPIHOSTFILE))

  def getConf: Configuration = {
    val conf: Configuration = confBroadcast.value
    log.info(conf.toString)
    log.info(conf.get(GenomicsDBConfiguration.LOADERJSON))
    if (shouldCloneJobConf) {
        /**
          * Picked up from [[org.apache.spark.rdd.NewHadoopRDD]]
          */
        GenomicsDBRDD.CONFIGURATION_INSTANTIATION_LOCK.synchronized {
        log.debug("Cloning Hadoop Configuration in GenomicsDBRDD")
        new Configuration(conf)
      }
    } else {
      conf.set(GenomicsDBConfiguration.LOADERJSON, loaderJsonBroadcast.value)
      conf.set(GenomicsDBConfiguration.QUERYJSON, queryJsonBroadcast.value)
      conf.set(GenomicsDBConfiguration.MPIHOSTFILE, hostsBroadcast.value)
      conf
    }
  }

  /**
    * Method containing the next() method to iterate over the variants. The first call
    * reads variants from underlying GenomicsDB into a internal buffer. Subsequent calls
    * read from the buffer until al l records are read. Then, next call to GenomicsDB is
    * executed.
    *
    * @param thisSplit  the local partition of the task. GenomicsDB inputsplit carries the
    *                   TileDB context
    * @param context    task context carries configuration
    * @return           returns the iterators over variants
    */
  @DeveloperApi
  @Override
  def compute(thisSplit: Partition, context: TaskContext): Iterator[VCONTEXT] = {

    val iterator = new Iterator[VCONTEXT] {
      val split: GenomicsDBPartition[VCONTEXT, SOURCE] =
        thisSplit.asInstanceOf[GenomicsDBPartition[VCONTEXT, SOURCE]]
      val conf: Configuration = getConf
      val inputFormat = new GenomicsDBInputFormat[VCONTEXT, SOURCE]
      inputFormat.setConf(conf)

      val recordReader: RecordReader[String, VCONTEXT] =
        inputFormat.createRecordReader(split.serializableSplit.value, null)
      val gRecordReader: GenomicsDBRecordReader[VCONTEXT, SOURCE] =
        recordReader.asInstanceOf[GenomicsDBRecordReader[VCONTEXT, SOURCE]]
      gRecordReader.initialize(null, null)

      var finished = false
      var havePair = false

      /**
        * Check whether more variants exist or not
        *
        * @return  true or false
        */
      override def hasNext: Boolean = {
        if (!finished && !havePair) {
          finished = !gRecordReader.nextKeyValue
          if (finished) {
            // Close and release the reader here; close() will also be called when the task
            // completes, but for tasks that read from many files, it helps to release the
            // resources early.
            close()
          }
          havePair = !finished
        }
        !finished
      }

      private def close() {
        gRecordReader.close()
      }

      override def next(): VCONTEXT = {
        if (!hasNext) {
          throw new NoSuchElementException("End of stream")
        }
        havePair = false
        gRecordReader.getCurrentValue
      }
    }
    new InterruptibleIterator[VCONTEXT](context, iterator)
  }

  /**
    * @todo  Job context is weird!! Make it go away
    * @return  return list of GenomicsDB partitions
    */
  @Override
  protected def getPartitions: Array[Partition] = {

    val inputFormat = new GenomicsDBInputFormat[VCONTEXT, SOURCE]
    inputFormat.setConf(sparkContext.hadoopConfiguration)
    val rawSplits = inputFormat.getSplits(null).toArray
    val result = new Array[Partition](rawSplits.length)
    for (i <- 0 until rawSplits.length) {
      result(i) = new GenomicsDBPartition(id, i,
        rawSplits(i).asInstanceOf[InputSplit with Writable])
    }
    result
  }

  /**
    * Return a collection of Java objects
    *
    * @return
    */
  override def toJavaRDD: JavaRDD[VCONTEXT] = super.toJavaRDD()

  override def persist(storageLevel: StorageLevel): this.type = {
    log.debug("GenomicsDB persist called with level: " + storageLevel.toString)
    if (storageLevel.equals(StorageLevel.DISK_ONLY) |
      storageLevel.equals(StorageLevel.MEMORY_AND_DISK) ) {
      super.persist
    }
    this
  }

  /**
    * Returns the count of the variants in GenomicsDB
    *
    * @return  integer count of variants
    */
  override def count(): Long = {
    super.count()
  }

  @Override override def collect(): Array[VCONTEXT] = {
    super.collect()
  }
}

object GenomicsDBRDD {
  val CONFIGURATION_INSTANTIATION_LOCK = new Object()
}