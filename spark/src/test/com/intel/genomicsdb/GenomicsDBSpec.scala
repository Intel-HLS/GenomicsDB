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

import org.scalatest.FunSuite

class GenomicsDBSpec extends FunSuite {

  test("Spark Standalone Example") {

    /* Create a GenomicsDB Context which will eventually create a Spark Context */
    val genomicsDBConfiguration = new GenomicsDBConf()
      .setMaster("spark://GAO-WS5:7077")
      .setAppName("GenomicsDB StandAlone Example")
      .setLoaderJsonFile("/home/kdatta1/GenomicsDB/GenomicsDB/tests/loader.json")
      .setQueryJsonFile("/home/kdatta1/GenomicsDB/GenomicsDB/tests/query.json")
      .setHostFile("/home/kdatta1/GenomicsDB/GenomicsDB/tests/hostfile")

    val genomicsDBContext = new GenomicsDBContext(genomicsDBConfiguration)
    val myRDD = genomicsDBContext.getVariantContexts
    myRDD.count()
  }

}
