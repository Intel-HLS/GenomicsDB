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

import htsjdk.tribble.Feature
import org.apache.hadoop.conf.Configuration
import org.apache.hadoop.io.Writable
import org.apache.hadoop.mapreduce.InputSplit
import org.apache.spark._
import org.apache.spark.annotation.DeveloperApi
import org.apache.spark.api.java.JavaRDD
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
  log.info("Just before calling confBroadcast::")
  private val confBroadcast = gc.sparkContext.broadcast(
    new GenomicsDBConf(gc.sparkContext.hadoopConfiguration))
  log.info("Just after calling confBroadcast::" + confBroadcast.value.toString)
  log.info(confBroadcast.value.get(GenomicsDBConf.LOADERJSON))
  log.info(confBroadcast.value.get(GenomicsDBConf.QUERYJSON))

  def getConf: Configuration = {
    val conf: Configuration = confBroadcast.value
    log.info(conf.toString)
    log.info(conf.get(GenomicsDBConf.LOADERJSON))
    if (shouldCloneJobConf) {
      // Picked up from NewHadoopRDD::
      // Hadoop Configuration objects are not thread-safe, which may lead to various problems if
      // one job modifies a configuration while another reads it (SPARK-2546, SPARK-10611).  This
      // problem occurs somewhat rarely because most jobs treat the configuration as though it's
      // immutable.  One solution, implemented here, is to clone the Configuration object.
      // Unfortunately, this clone can be very expensive.  To avoid unexpected performance
      // regressions for workloads and Hadoop versions that do not suffer from these thread-safety
      // issues, this cloning is disabled by default.
      GenomicsDBRDD.CONFIGURATION_INSTANTIATION_LOCK.synchronized {
        log.debug("Cloning Hadoop Configuration in GenomicsDBRDD")
        // The Configuration passed in is actually a JobConf and possibly contains credentials.
        // To keep those credentials properly we have to create a new JobConf not a Configuration.

        log.info("getConf::77::" + conf.get(GenomicsDBConf.LOADERJSON))
        log.info("getConf::78::" + conf.get(GenomicsDBConf.QUERYJSON))
        new Configuration(conf)
      }
    } else {
      log.info("getConf::83::" + conf.get(GenomicsDBConf.LOADERJSON))
      log.info("getConf::84::" + conf.get(GenomicsDBConf.QUERYJSON))
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

    log.info("compute() called")
    val iterator = new Iterator[VCONTEXT] {
      val split = thisSplit.asInstanceOf[GenomicsDBPartition[VCONTEXT, SOURCE]]
      val conf = getConf
      val inputFormat = new GenomicsDBInputFormat[VCONTEXT, SOURCE]
      inputFormat.setConf(conf)

      val recordReader = inputFormat.createRecordReader(split.serializableSplit.value, null)
      val gRecordReader =
        recordReader.asInstanceOf[GenomicsDBRecordReader[VCONTEXT, SOURCE]]
      gRecordReader.initialize(null, null)

      var havePair = false
      var finished = false

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
        if (gRecordReader!=null) {
          // Note: it's important that we set closed = true before calling close(), since setting it
          // afterwards would permit us to call close() multiple times if close() threw an exception.
          gRecordReader.close()
        }
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
    log.info("getPartitions called")
    val inputFormat = new GenomicsDBInputFormat[VCONTEXT, SOURCE]
    inputFormat.setConf(confBroadcast.value)
    val rawSplits = inputFormat.getSplits(null).toArray
    val result = new Array[Partition](rawSplits.length)
    for (i <- 0 until rawSplits.length) {
      result(i) = new GenomicsDBPartition(id, i,
        rawSplits(i).asInstanceOf[InputSplit with Writable])
    }
    log.info("getPartitions done")
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
  override def count() = {
    super.count()
  }

  @Override override def collect(): Array[VCONTEXT] = {
    super.collect()
  }
}

object GenomicsDBRDD {
  val CONFIGURATION_INSTANTIATION_LOCK = new Object()
}