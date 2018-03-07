/*
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

package com.intel.genomicsdb;

import com.intel.genomicsdb.model.GenomicsDBCallsetsMapProto;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;

public final class GenomicsDBCallsetsMapProtoSpec {

  @Test(testName = "Validate callset map protocol buffer")
  public void testCallsetMapProto() {
    GenomicsDBCallsetsMapProto.CallsetMappingPB.Builder cBuilder =
      GenomicsDBCallsetsMapProto.CallsetMappingPB.newBuilder();

    GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap.Builder sampleBuilder0 =
      GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap.newBuilder();

    GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap sample0 =
      sampleBuilder0
        .setSampleName("ABC")
        .setIdxInFile(0)
        .setRowIdx(0)
        .build();

    GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap.Builder sampleBuilder1 =
      GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap.newBuilder();

    GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap sample1 =
      sampleBuilder1
        .setSampleName("MNP")
        .setIdxInFile(0)
        .setRowIdx(1)
        .build();

    GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap.Builder sampleBuilder2 =
      GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap.newBuilder();

    GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap sample2 =
      sampleBuilder2
        .setSampleName("XYZ")
        .setIdxInFile(0)
        .setRowIdx(2)
        .build();

    GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB =
      cBuilder
        .addCallsets(sample0)
        .addCallsets(sample1)
        .addCallsets(sample2)
        .build();


    Assert.assertEquals(callsetMappingPB.isInitialized(), true);
    List<GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap> sampleList =
      callsetMappingPB.getCallsetsList();

    Assert.assertEquals(sample0, sampleList.get(0));
    Assert.assertEquals(sample1, sampleList.get(1));
    Assert.assertEquals(sample2, sampleList.get(2));
  }
}
