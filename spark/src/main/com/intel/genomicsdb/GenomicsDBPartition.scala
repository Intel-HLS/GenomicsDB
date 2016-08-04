package com.intel.genomicsdb

import htsjdk.tribble.Feature
import org.apache.hadoop.mapred.InputSplit
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
    rawSplit: InputSplit)
  extends Partition {

  val serializableSplit = new SerializableWritable[InputSplit](rawSplit)

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