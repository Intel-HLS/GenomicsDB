package com.intel.genomicsdb

import genomicsdb.GenomicsDBFeatureReader
import htsjdk.tribble.{CloseableTribbleIterator, Feature}
import org.apache.hadoop.io.LongWritable
import org.apache.hadoop.mapreduce.{InputSplit, RecordReader, TaskAttemptContext}

import scala.reflect.ClassTag

/**
  * A RecordReader class to read variant records from GenomicsDB.
  * It emits serializable vatiant contexts so that RDDs can use
  * them. The type of keys are NULL, and the type of values are
  * VariantContexts
  *
  * @param split
  * @param conf
  * @tparam D
  * @tparam S
  */
class GenomicsDBRecordReader[D <: Feature: ClassTag, S: ClassTag](
    split: GenomicsDBInputSplit,
    conf: GenomicsDBConf,
    featureReader: GenomicsDBFeatureReader[D, S]) extends RecordReader[LongWritable, D] {

  var currentVariant : D = _
  var iterator: CloseableTribbleIterator[D] = _

  @Override
  def getProgress: Float = {
    0
  }

  //def nextKeyValue: D = iterator.next()

  @Override
  def getCurrentValue: D = {
    currentVariant
  }

  @Override
  def initialize(inputSplit: InputSplit, taskAttemptContext: TaskAttemptContext): Unit = {
    initialize
  }

  def initialize: Unit = {
    iterator = featureReader.iterator
    if (hasNext)
      currentVariant = iterator.next
  }

  @Override
  def getCurrentKey: LongWritable = new LongWritable(currentVariant.toString.hashCode.toLong)

  @Override
  def close(): Unit = {
    iterator.close()
  }

  override def nextKeyValue(): Boolean = {
    if (hasNext) true else false
  }

  def nextValue: D = iterator.next

  def hasNext: Boolean = iterator.hasNext
}
