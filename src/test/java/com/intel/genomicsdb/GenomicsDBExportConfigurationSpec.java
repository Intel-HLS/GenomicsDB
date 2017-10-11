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

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GenomicsDBExportConfigurationSpec {

  private static final String TEST_ARRAY = "array";
  private final String TILEDB_WORKSPACE = "/path/to/junk/folder";
  private final String REFERENCE_GENOME = "/path/to/HG38.fasta";
  private final String TEST_VIDMAP_FILENAME = "/path/to/junk/vidmap.json";
  private final String TEST_CALLSETMAP_FILENAME = "/path/to/junk/callsetmap.json";
  private final List<String> ATTRIBUTES = new ArrayList<>(Arrays.asList("REF", "ALT", "BaseQRankSum", "MQ", "RAW_MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP", "PID", "PGT"));

  @Test(groups = {"configuration tests"})
  public void testExportConfiguration() {
    List<GenomicsDBExportConfiguration.ColumnRange> columnRanges0 = new ArrayList<>(2);
    List<GenomicsDBExportConfiguration.ColumnRange> columnRanges1 = new ArrayList<>(2);
    List<GenomicsDBExportConfiguration.ColumnRangeList> columnRangeLists = new ArrayList<>(2);

    GenomicsDBExportConfiguration.ColumnRange.Builder columnRange0 =
      GenomicsDBExportConfiguration.ColumnRange.newBuilder();
    GenomicsDBExportConfiguration.ColumnRange r0 =
      columnRange0
        .setLow(0)
        .setHigh(100000)
        .build();
    
    GenomicsDBExportConfiguration.ColumnRange.Builder columnRange1 =
      GenomicsDBExportConfiguration.ColumnRange.newBuilder();
    GenomicsDBExportConfiguration.ColumnRange r1 =
      columnRange1
        .setLow(100001)
        .setHigh(200000)
        .build();

    columnRanges0.add(r0);
    columnRanges0.add(r1); 

    GenomicsDBExportConfiguration.ColumnRange.Builder columnRange2 =
      GenomicsDBExportConfiguration.ColumnRange.newBuilder();
    GenomicsDBExportConfiguration.ColumnRange r2 =
      columnRange2
        .setLow(200001)
        .setHigh(300000)
        .build();
    
    GenomicsDBExportConfiguration.ColumnRange.Builder columnRange3 =
      GenomicsDBExportConfiguration.ColumnRange.newBuilder();
    GenomicsDBExportConfiguration.ColumnRange r3 =
      columnRange3
        .setLow(300001)
        .setHigh(400000)
        .build();

    columnRanges1.add(r2);
    columnRanges1.add(r3); 

    GenomicsDBExportConfiguration.ColumnRangeList.Builder columnRangeList0 =
      GenomicsDBExportConfiguration.ColumnRangeList.newBuilder();

    GenomicsDBExportConfiguration.ColumnRangeList queryColumnRanges0 =
      columnRangeList0
        .addAllRangeList(columnRanges0)
        .build();

    GenomicsDBExportConfiguration.ColumnRangeList.Builder columnRangeList1 =
      GenomicsDBExportConfiguration.ColumnRangeList.newBuilder();

    GenomicsDBExportConfiguration.ColumnRangeList queryColumnRanges1 =
      columnRangeList1
        .addAllRangeList(columnRanges1)
        .build();

    columnRangeLists.add(queryColumnRanges0);
    columnRangeLists.add(queryColumnRanges1);


    GenomicsDBExportConfiguration.ExportConfiguration.Builder configBuilder =
      GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();
    GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration =
      configBuilder
        .setWorkspace(TILEDB_WORKSPACE)
        .setArray(TEST_ARRAY)
	.setReferenceGenome(REFERENCE_GENOME)
        .addAllQueryColumnRanges(columnRangeLists)        
	.addAllAttributes(ATTRIBUTES)
	.setVidMappingFile(TEST_VIDMAP_FILENAME)
	.setCallsetMappingFile(TEST_CALLSETMAP_FILENAME)
        .build();
    
    Assert.assertEquals(exportConfiguration.isInitialized(), true);

    // Assert has methods
    Assert.assertEquals(exportConfiguration.hasWorkspace(), true);
    Assert.assertEquals(exportConfiguration.hasArray(), true);
    Assert.assertEquals(exportConfiguration.hasReferenceGenome(), true);
    Assert.assertEquals(exportConfiguration.hasVidMappingFile(), true);
    Assert.assertEquals(exportConfiguration.hasCallsetMappingFile(), true);
    
    // Assert gets
    Assert.assertEquals(exportConfiguration.getWorkspace(), "/path/to/junk/folder");
    Assert.assertEquals(exportConfiguration.getArray(), "array");
    Assert.assertEquals(exportConfiguration.getReferenceGenome(), "/path/to/HG38.fasta");
    Assert.assertEquals(exportConfiguration.getQueryColumnRangesCount(), 2);
    Assert.assertEquals(exportConfiguration.getAttributesCount(), 19);
    Assert.assertEquals(exportConfiguration.getVidMappingFile(), TEST_VIDMAP_FILENAME);
    Assert.assertEquals(exportConfiguration.getCallsetMappingFile(), TEST_CALLSETMAP_FILENAME);
  }
}
