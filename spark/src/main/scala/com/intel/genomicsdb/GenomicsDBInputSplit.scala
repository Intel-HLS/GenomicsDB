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

import java.io.{DataInput, DataOutput}

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFHeader
import org.apache.commons.logging.{Log, LogFactory}
import org.apache.hadoop.classification.{InterfaceAudience, InterfaceStability}
import org.apache.hadoop.io.Writable
import org.apache.hadoop.mapreduce.InputSplit

/**
  * An InputSplit that spans a set of rows/cols in TileDB as defined by the cell ordering
  */
@InterfaceAudience.Public
@InterfaceStability.Evolving
class GenomicsDBInputSplit
  extends InputSplit with Writable {

  val log: Log = LogFactory.getLog("GenomicsDBInputSplit")
  private var hosts: Array[String] = _
  private var gConf: GenomicsDBConf = _

  /**
    * Return the locations of this split.
    * Although in Hadoop DFS, there can be multiple splits
    * of a block, in GenomicsDB, there's only one, for now!
    *
    * @return  Return the location as a singleton array
    */
  @Override
  def getLocations: Array[String] = {
    hosts
  }

  @Override
  def getLength: Long = 0

  /**
    * @todo  Implement a GenomicsDB handler which contains the TileDB storage manager
    * instance. Then, TileDB arrays can be directly written to
    *
    * @param header  VCF Header required to write data back to GenomicsDB
    * @param variants  VCF entry to be writen back to GenomicsDB
    */
  @Override
  def write(dataOutput: DataOutput, header: VCFHeader, variants: Array[VariantContext]): Unit = {
    dataOutput.writeChars(header.toString)
    for (i <- 0 to variants.length) {
      dataOutput.writeChars(variants(i).toString)
    }
  }

  @Override
  def readFields(dataInput: DataInput): Unit = {

  }

  @Override
  def write(dataOutput: DataOutput): Unit = {}

}

object GenomicsDBInputSplit {
  def apply(_hosts: Array[String], gconf: GenomicsDBConf): GenomicsDBInputSplit = {
    val split = new GenomicsDBInputSplit
    split.hosts = _hosts.clone()
    split.gConf = gconf.clone
    split
  }
}