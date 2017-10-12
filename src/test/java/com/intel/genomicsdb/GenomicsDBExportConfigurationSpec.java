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

  private final String VCF_HEADER_FILENAME = "inputs/template_vcf_header.vcf";
  private final String VCF_OUTPUT_FILENAME = "outputs/test.vcf.gz";
  private final String VCF_OUTPUT_FORMAT = "z0";
  private final int MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED = 50;
  private final boolean INDEX_OUTPUT_VCF = true;

  @Test(groups = {"configuration tests"})
  public void testExportConfiguration() {
    List<GenomicsDBExportConfiguration.ColumnRange> columnRanges0 = new ArrayList<>(2);
    List<GenomicsDBExportConfiguration.ColumnRange> columnRanges1 = new ArrayList<>(2);
    List<GenomicsDBExportConfiguration.ColumnRangeList> columnRangeLists = new ArrayList<>(2);
    List<GenomicsDBExportConfiguration.RowRange> rowRanges0 = new ArrayList<>(2);
    List<GenomicsDBExportConfiguration.RowRange> rowRanges1 = new ArrayList<>(2);
    List<GenomicsDBExportConfiguration.RowRangeList> rowRangeLists = new ArrayList<>(2);


    GenomicsDBExportConfiguration.ColumnRange.Builder columnRange0 =
      GenomicsDBExportConfiguration.ColumnRange.newBuilder();
    GenomicsDBExportConfiguration.ColumnRange c0 =
      columnRange0
        .setLow(0)
        .setHigh(100000)
        .build();
    
    GenomicsDBExportConfiguration.ColumnRange.Builder columnRange1 =
      GenomicsDBExportConfiguration.ColumnRange.newBuilder();
    GenomicsDBExportConfiguration.ColumnRange c1 =
      columnRange1
        .setLow(100001)
        .setHigh(200000)
        .build();

    columnRanges0.add(c0);
    columnRanges0.add(c1); 

    GenomicsDBExportConfiguration.ColumnRange.Builder columnRange2 =
      GenomicsDBExportConfiguration.ColumnRange.newBuilder();
    GenomicsDBExportConfiguration.ColumnRange c2 =
      columnRange2
        .setLow(200001)
        .setHigh(300000)
        .build();
    
    GenomicsDBExportConfiguration.ColumnRange.Builder columnRange3 =
      GenomicsDBExportConfiguration.ColumnRange.newBuilder();
    GenomicsDBExportConfiguration.ColumnRange c3 =
      columnRange3
        .setLow(300001)
        .setHigh(400000)
        .build();

    columnRanges1.add(c2);
    columnRanges1.add(c3); 

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

    GenomicsDBExportConfiguration.RowRange.Builder rowRange0 =
      GenomicsDBExportConfiguration.RowRange.newBuilder();
    GenomicsDBExportConfiguration.RowRange r0 =
      rowRange0
        .setLow(0)
        .setHigh(100)
        .build();
    
    GenomicsDBExportConfiguration.RowRange.Builder rowRange1 =
      GenomicsDBExportConfiguration.RowRange.newBuilder();
    GenomicsDBExportConfiguration.RowRange r1 =
      rowRange1
        .setLow(101)
        .setHigh(200)
        .build();

    rowRanges0.add(r0);
    rowRanges0.add(r1); 

    GenomicsDBExportConfiguration.RowRange.Builder rowRange2 =
      GenomicsDBExportConfiguration.RowRange.newBuilder();
    GenomicsDBExportConfiguration.RowRange r2 =
      rowRange2
        .setLow(201)
        .setHigh(300)
        .build();
    
    GenomicsDBExportConfiguration.RowRange.Builder rowRange3 =
      GenomicsDBExportConfiguration.RowRange.newBuilder();
    GenomicsDBExportConfiguration.RowRange r3 =
      rowRange3
        .setLow(301)
        .setHigh(400)
        .build();

    rowRanges1.add(r2);
    rowRanges1.add(r3); 

    GenomicsDBExportConfiguration.RowRangeList.Builder rowRangeList0 =
      GenomicsDBExportConfiguration.RowRangeList.newBuilder();

    GenomicsDBExportConfiguration.RowRangeList queryRowRanges0 =
      rowRangeList0
        .addAllRangeList(rowRanges0)
        .build();

    GenomicsDBExportConfiguration.RowRangeList.Builder rowRangeList1 =
      GenomicsDBExportConfiguration.RowRangeList.newBuilder();

    GenomicsDBExportConfiguration.RowRangeList queryRowRanges1 =
      rowRangeList1
        .addAllRangeList(rowRanges1)
        .build();

    rowRangeLists.add(queryRowRanges0);
    rowRangeLists.add(queryRowRanges1);


    GenomicsDBExportConfiguration.ExportConfiguration.Builder configBuilder =
      GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();
    GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration =
      configBuilder
        .setWorkspace(TILEDB_WORKSPACE)
        .setArray(TEST_ARRAY)
	.setReferenceGenome(REFERENCE_GENOME)
        .addAllQueryColumnRanges(columnRangeLists)        
        .addAllQueryRowRanges(rowRangeLists)        
	.addAllAttributes(ATTRIBUTES)
	.setVidMappingFile(TEST_VIDMAP_FILENAME)
	.setCallsetMappingFile(TEST_CALLSETMAP_FILENAME)
	.setVcfHeaderFilename(VCF_HEADER_FILENAME)
	.setVcfOutputFilename(VCF_OUTPUT_FILENAME)
	.setVcfOutputFormat(VCF_OUTPUT_FORMAT)
	.setMaxDiploidAltAllelesThatCanBeGenotyped(MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED)
	.setIndexOutputVCF(INDEX_OUTPUT_VCF)
        .build();
    
    Assert.assertEquals(exportConfiguration.isInitialized(), true);

    // Assert has methods
    Assert.assertEquals(exportConfiguration.hasWorkspace(), true);
    Assert.assertEquals(exportConfiguration.hasArray(), true);
    Assert.assertEquals(exportConfiguration.hasReferenceGenome(), true);
    Assert.assertEquals(exportConfiguration.hasVidMappingFile(), true);
    Assert.assertEquals(exportConfiguration.hasCallsetMappingFile(), true);
    Assert.assertEquals(exportConfiguration.hasVcfHeaderFilename(), true);
    Assert.assertEquals(exportConfiguration.hasVcfOutputFilename(), true);                                              
    Assert.assertEquals(exportConfiguration.hasVcfOutputFormat(), true);
    Assert.assertEquals(exportConfiguration.hasMaxDiploidAltAllelesThatCanBeGenotyped(), true);
    Assert.assertEquals(exportConfiguration.hasIndexOutputVCF(), true);
    
    // Assert gets
    Assert.assertEquals(exportConfiguration.getWorkspace(), TILEDB_WORKSPACE);
    Assert.assertEquals(exportConfiguration.getArray(), TEST_ARRAY);
    Assert.assertEquals(exportConfiguration.getReferenceGenome(), REFERENCE_GENOME);
    Assert.assertEquals(exportConfiguration.getQueryColumnRangesCount(), 2);
    Assert.assertEquals(exportConfiguration.getQueryRowRangesCount(), 2);
    Assert.assertEquals(exportConfiguration.getAttributesCount(), 19);
    Assert.assertEquals(exportConfiguration.getVidMappingFile(), TEST_VIDMAP_FILENAME);
    Assert.assertEquals(exportConfiguration.getCallsetMappingFile(), TEST_CALLSETMAP_FILENAME);
    Assert.assertEquals(exportConfiguration.getVcfHeaderFilename(),VCF_HEADER_FILENAME);
    Assert.assertEquals(exportConfiguration.getVcfOutputFilename(),VCF_OUTPUT_FILENAME);                                              
    Assert.assertEquals(exportConfiguration.getVcfOutputFormat(),VCF_OUTPUT_FORMAT);
    Assert.assertEquals(exportConfiguration.getMaxDiploidAltAllelesThatCanBeGenotyped(),MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED);
    Assert.assertEquals(exportConfiguration.getIndexOutputVCF(),INDEX_OUTPUT_VCF);
  }
}
