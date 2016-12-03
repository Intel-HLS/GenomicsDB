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

import java.io.{DataOutput, FileInputStream}
import java.util
import java.util.Scanner

import org.apache.hadoop.classification.{InterfaceAudience, InterfaceStability}
import org.apache.hadoop.mapreduce.{Job, JobContext}
import org.apache.spark.SparkConf

import scala.reflect.io.Path

@InterfaceAudience.Public
@InterfaceStability.Evolving
class GenomicsDBConf extends SparkConf with Serializable {

  import GenomicsDBConf._

  val hosts = new util.ArrayList[String]

  /**
    * Set the minimum input split size
    *
    * @param job  the job to modify
    * @param size the minimum size
    */
  def setMinInputSplitSize(job: Job, size: Long) {
    job.getConfiguration.setLong(SPLIT_MINSIZE, size)
  }

  /**
    * Get the minimum split size
    *
    * @param job the job
    * @return the minimum number of bytes that can be in a split
    */
  def getMinSplitSize(job: JobContext): Long = job.getConfiguration.getLong(SPLIT_MINSIZE, 1L)

  /**
    * Set the 64-bit value for a configuration parameter
    *
    * @param key
    * @param value
    * @return
    */
  def setLong(key: String, value: Long): GenomicsDBConf = {
    set(key, value.toString)
    this
  }

  /**
    * Return the 64-bit integer value
    *
    * @param key
    * @return
    */
  def getLong(key: String): Long = get(key).toLong


  /**
    * Set the maximum split size
    *
    * @param job  the job to modify
    * @param size the maximum split size
    */
  def setMaxInputSplitSize(job: Job, size: Long): GenomicsDBConf = {
    setLong(SPLIT_MAXSIZE, size)
    this
  }

  /**
    * Get the maximum split size.
    *
    * @param context the job to look at.
    * @return the maximum number of bytes a split can include
    */
  def getMaxSplitSize(context: JobContext): Long =
    getLong(SPLIT_MAXSIZE, Long.MaxValue)

  /**
    * Set the loader JSON file if one is provided
    *
    * @param loaderFile  path of the JSON loader file
    */
  def setLoaderJsonFile(loaderFile: Path): GenomicsDBConf = {
    set(LOADERJSON, loaderFile.toString)
    this
  }

  def getLoaderJson: String = get(LOADERJSON)

  /**
    * Set the query JSON file if one is provided
    *
    * @param queryFile  path of the JSON query file
    */
  def setQueryJsonFile(queryFile: Path): GenomicsDBConf = {
    set(QUERYJSON, queryFile.toString)
    this
  }

  def getQueryJson: String = get(QUERYJSON)

  /**
    * Set the host file if one is provided
    *
    * @param file  full path of the MPI host file
    */
  def setHostFile(file: Path): GenomicsDBConf = {
    set(MPIHOSTFILE, file.toString)
    val scanner = new Scanner(new FileInputStream(file.toString))
    while (scanner.hasNextLine) {
      hosts.add(scanner.nextLine())
    }
    scanner.close()
    this
  }

  def getHostFile = get(MPIHOSTFILE)

  override def setMaster(value: String): GenomicsDBConf = {
    super.setMaster(value)
    this
  }

  override def setAppName(value: String): GenomicsDBConf = {
    super.setAppName(value)
    this
  }

  def getHosts: util.ArrayList[String] = hosts

  def setGenomicsDBPort(port: Long): GenomicsDBConf = {
    setLong(THRIFT_PORT, port)
    this
  }

  def getGenomicsDBPort: Int = getLong(THRIFT_PORT).toInt

  def setTileDBWorkspace(directory: Path): GenomicsDBConf = {
    set(TILEDB_WORKSPACE, directory.toString)
    this
  }

  def getTileDBWorkspace: Path = get(TILEDB_WORKSPACE)

  def setTileDBArrayName(name: String) : GenomicsDBConf = {
    set(TILEDB_ARRAYNAME, name)
    this
  }

  def getTileDBArrayName: String = get(TILEDB_ARRAYNAME)

  def setReferenceGenome(file: String): GenomicsDBConf = {
    set(REFERENCE_GENOME_FILE, file)
    this
  }

  def getReferenceGenome: String = get(REFERENCE_GENOME_FILE)

  def setVCFHeaderFile(file: String): GenomicsDBConf = {
    set(VCF_HEADER_FILE, file)
    this
  }

  def getVCFHeaderFile: String = get(VCF_HEADER_FILE)

  @Override
  def write(dataOutput: DataOutput): Unit = {
    for (i <- 0 until hosts.size())
      dataOutput.writeChars(hosts.get(i)+"\n")
  }

  override def toString(): String = {
    val str = new StringBuilder
    for (i <- 0 until hosts.size())
      str.+(hosts.get(i) + "\n")
    str.toString()
  }

  override def clone(): GenomicsDBConf = {
    this
  }

//  @Override
//  def readFields(dataInput: DataInput): Unit = {
//    for (i <- hosts.size()) {
//      hosts.add(dataInput.readLine())
//    }
//  }
}


object GenomicsDBConf {

  val TILEDB_WORKSPACE: String = "genomicsdb.inputformat.workspace"
  val TILEDB_ARRAYNAME: String = "genomicsdb.input.array.name"
  val SPLIT_MINSIZE: String = "genomicsdb.inputformat.split.minsize"
  val SPLIT_MAXSIZE: String = "genomicsdb.inputformat.split.maxsize"
  val LOADERJSON: String = "genomicsdb.input.loaderjsonfile"
  val QUERYJSON: String = "genomisdb.input.queryjsonfile"
  val MPIHOSTFILE: String = "genomisdb.input.mpi.hostfile"
  val THRIFT_PORT: String = "genomicsdb.input.thrift.port"
  val REFERENCE_GENOME_FILE: String = "genomicsdb.input.reference.genome"
  val VCF_HEADER_FILE: String = "genomicsdb.input.vcf.header.file"
}
