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

import htsjdk.tribble.readers.PositionalBufferedStream
import htsjdk.variant.variantcontext.VariantContext
import org.apache.hadoop.conf.Configuration
import org.apache.log4j.Logger
import org.apache.spark.{SparkConf, SparkContext}
import org.apache.spark.rdd.RDD

/**
  * GenomicsDB context
  *
  * @param conf  Hadoop configuration including GenomicsDB
  *              keys. Copy these keys in case we have to
  *              create a new SparkContext
  */
class GenomicsDBContext(val conf: Configuration, var sparkContext: SparkContext = null) {

  val logger = Logger.getLogger(classOf[GenomicsDBContext])

  if (sparkContext == null) {
    val conf = new SparkConf()
    conf.setMaster(conf.get("spark.master"))
    conf.setAppName("GenomicsDBTest")

    sparkContext = new SparkContext(conf)
    logger.info("new Spark Context created from GenomicsDB Context")
  } else {
    logger.info("Creating GenomicsDB Context from Spark Context")
  }

  def getVariantContexts: RDD[VariantContext] = {
    new GenomicsDBRDD[VariantContext, PositionalBufferedStream](this, 1)
  }
}

