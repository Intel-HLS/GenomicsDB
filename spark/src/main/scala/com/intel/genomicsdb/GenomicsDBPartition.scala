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
import org.apache.hadoop.io.Writable
import org.apache.hadoop.mapreduce.InputSplit
import org.apache.spark.{Partition, SerializableWritable}

/**
  * GenomicsDBPartition is a container class for RDD partitions for the underlying variant
  * data. It has a 1:1 mapping between GenomicsDB input splits
  *
  * @param rddId  unique rdd identifier
  * @param index  index for task scheduler
  * @param rawSplit  GenomicsDB input split passed from RDD getpartitions
  *                  [[com.intel.genomicsdb.GenomicsDBInputSplit]]
  */
private[genomicsdb] class GenomicsDBPartition[VCONTEXT <: Feature, SOURCE](
    rddId: Int,
    val index: Int,
    rawSplit: InputSplit with Writable)
  extends Partition {

  val serializableSplit = new SerializableWritable(rawSplit)

  override def hashCode: Int = 31 * (31 + rddId) + index

  override def equals(other: Any): Boolean = super.equals(other)
}

/**
  * Partitioning strategy can be either row-major or column-major
  * The value is specified in partitioning schema in loader JSON
  * file
  */
object GenomicsDBPartitionStrategy extends Enumeration {
  val ROW_MAJOR, COL_MAJOR = Value
}