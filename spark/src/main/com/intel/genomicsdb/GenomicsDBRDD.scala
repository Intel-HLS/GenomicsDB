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
import org.apache.hadoop.io.{LongWritable, Writable}
import org.apache.hadoop.mapreduce.{InputSplit, RecordReader}
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

  val gConf = gc.gconf

  // Loader and query configuration scripts are meant to be about few megabytes,
  // be careful about this although it still makes sense to broadcast it as
  // they are used in C++ GenomicsDB interface
  private val loaderBroadcast = gc.sparkContext.broadcast(gConf.getLoaderJson)
  private val queryBroadcast = gc.sparkContext.broadcast(gConf.getQueryJson)
  private val confBroadcast = gc.sparkContext.broadcast(gConf)

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
    val split = thisSplit.asInstanceOf[GenomicsDBPartition[VCONTEXT, SOURCE]]
    log.debug("Inside GenomicsDB.compute")
    System.out.println("Inside GenomicsDB.compute")
    val inputFormat = new GenomicsDBInputFormat[VCONTEXT, SOURCE](gConf)
    inputFormat.setLoaderJsonFile(gConf.getLoaderJson)
    inputFormat.setQueryJsonFile(gConf.getQueryJson)
    inputFormat.setHostFile(gConf.getHostFile)

    val iterator = new Iterator[VCONTEXT] {

      // Register an on-task-completion callback to close the input stream.
      //context.addTaskCompletionListener{ context => closeIfNeeded() }

      private var gotNext = false
      private var nextValue: VCONTEXT = _
      private var closed = false
      protected var finished = false

      val recordReader: RecordReader[LongWritable, VCONTEXT] =
        inputFormat.createRecordReader(split.serializableSplit.value,
          context, 0)

      override def hasNext: Boolean = {
        if (!finished) {
          if (!gotNext) {
            nextValue = getNext
            if (finished) {
              closeIfNeeded()
            }
            gotNext = true
          }
        }
        !finished
      }

      def closeIfNeeded() {
        if (!closed) {
          // Note: it's important that we set closed = true before calling close(), since setting it
          // afterwards would permit us to call close() multiple times if close() threw an exception.
          closed = true
          recordReader.close()
        }
      }

      override def next(): VCONTEXT = {
        if (!hasNext) {
          throw new NoSuchElementException("End of stream")
        }
        gotNext = false
        nextValue
      }

      def getNext: VCONTEXT = {
        recordReader.nextKeyValue()
        recordReader.getCurrentValue
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
    log.debug("Getting partitions for GenomicsDBRDD")
    System.out.println("Getting partitions for GenomicsDBRDD")
    val jobConf = confBroadcast
    val inputFormat = new GenomicsDBInputFormat[VCONTEXT, SOURCE](gc.gconf)
    val rawSplits = inputFormat.getSplits(runningAsService)
    val result = new Array[Partition](rawSplits.size)
    for (i <- rawSplits.indices) {
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
}