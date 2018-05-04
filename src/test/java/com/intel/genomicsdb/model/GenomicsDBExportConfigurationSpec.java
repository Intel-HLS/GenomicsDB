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
package com.intel.genomicsdb.model;

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
  private final List<String> ATTRIBUTES = new ArrayList<>(Arrays.asList("REF", "ALT", "BaseQRankSum", "MQ",
          "RAW_MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL",
          "DP_FORMAT", "MIN_DP", "PID", "PGT"));

  private final String VCF_HEADER_FILENAME = "inputs/template_vcf_header.vcf";
  private final String VCF_OUTPUT_FILENAME = "outputs/test.vcf.gz";
  private final String VCF_OUTPUT_FORMAT = "z0";
  private final int MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED = 50;
  private final boolean INDEX_OUTPUT_VCF = true;

  @Test(groups = {"configuration tests"})
  public void testExportConfiguration() {
    List<GenomicsDBExportConfiguration.RowRange> rowRanges0 = new ArrayList<>(2);
    List<GenomicsDBExportConfiguration.RowRange> rowRanges1 = new ArrayList<>(2);
    List<GenomicsDBExportConfiguration.RowRangeList> rowRangeLists = new ArrayList<>(2);

    Coordinates.ContigInterval interval0 = Coordinates.ContigInterval.newBuilder().setBegin(0).setEnd(100000).setContig("").build();
    Coordinates.ContigInterval interval1 = Coordinates.ContigInterval.newBuilder().setBegin(100001).setEnd(200000).setContig("").build();
    Coordinates.GenomicsDBColumnInterval columnInterval0 = Coordinates.GenomicsDBColumnInterval.newBuilder()
            .setContigInterval(interval0).setContigInterval(interval1).build();
    Coordinates.GenomicsDBColumnOrInterval columnRanges0 = Coordinates.GenomicsDBColumnOrInterval.newBuilder()
            .setColumnInterval(columnInterval0).build();

    Coordinates.ContigInterval interval2 = Coordinates.ContigInterval.newBuilder().setBegin(200001).setEnd(300000).setContig("").build();
    Coordinates.GenomicsDBColumnInterval columnInterval1 = Coordinates.GenomicsDBColumnInterval.newBuilder()
            .setContigInterval(interval1).setContigInterval(interval2).build();
    Coordinates.GenomicsDBColumnOrInterval columnRanges1 = Coordinates.GenomicsDBColumnOrInterval.newBuilder()
            .setColumnInterval(columnInterval1).build();

  List<Coordinates.GenomicsDBColumnOrInterval> columnRangeLists = new ArrayList<>(2);
  columnRangeLists.add(columnRanges0);
  columnRangeLists.add(columnRanges1);

  GenomicsDBExportConfiguration.GenomicsDBColumnOrIntervalList list =
          GenomicsDBExportConfiguration.GenomicsDBColumnOrIntervalList.newBuilder()
                  .addAllColumnOrIntervalList(columnRangeLists).build();


    GenomicsDBExportConfiguration.RowRange.Builder rowRange0 = GenomicsDBExportConfiguration.RowRange.newBuilder();
    GenomicsDBExportConfiguration.RowRange r0 = rowRange0.setLow(0).setHigh(100).build();
    
    GenomicsDBExportConfiguration.RowRange.Builder rowRange1 = GenomicsDBExportConfiguration.RowRange.newBuilder();
    GenomicsDBExportConfiguration.RowRange r1 = rowRange1.setLow(101).setHigh(200).build();

    rowRanges0.add(r0);
    rowRanges0.add(r1); 

    GenomicsDBExportConfiguration.RowRange.Builder rowRange2 = GenomicsDBExportConfiguration.RowRange.newBuilder();
    GenomicsDBExportConfiguration.RowRange r2 = rowRange2.setLow(201).setHigh(300).build();
    
    GenomicsDBExportConfiguration.RowRange.Builder rowRange3 = GenomicsDBExportConfiguration.RowRange.newBuilder();
    GenomicsDBExportConfiguration.RowRange r3 = rowRange3.setLow(301).setHigh(400).build();

    rowRanges1.add(r2);
    rowRanges1.add(r3); 

    GenomicsDBExportConfiguration.RowRangeList.Builder rowRangeList0 = GenomicsDBExportConfiguration.RowRangeList.newBuilder();

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
        .setArrayName(TEST_ARRAY)
    .setReferenceGenome(REFERENCE_GENOME)
        .addQueryColumnRanges(list)
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
    Assert.assertEquals(exportConfiguration.hasArrayName(), true);
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
    Assert.assertEquals(exportConfiguration.getArrayName(), TEST_ARRAY);
    Assert.assertEquals(exportConfiguration.getReferenceGenome(), REFERENCE_GENOME);
    Assert.assertEquals(exportConfiguration.getQueryColumnRangesCount(), 1);

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
