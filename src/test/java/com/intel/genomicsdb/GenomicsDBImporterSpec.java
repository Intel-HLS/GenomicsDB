/*
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

package com.intel.genomicsdb;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class GenomicsDBImporterSpec {

  private static final String TILEDB_WORKSPACE = "./__workspace";
  private static final String ARRAY_FOR_PARTITION0 = "array_test0";
  private static final String ARRAY_FOR_PARTITION1 = "array_test1";

  @Test(groups = {"genomicsdb importer with an interval and multiple GVCFs"})
  public void testMultiGVCFInputs() throws IOException {

    File t6 = new File("tests/inputs/vcfs/t6.vcf.gz");
    File t7 = new File("tests/inputs/vcfs/t7.vcf.gz");
    File t8 = new File("tests/inputs/vcfs/t8.vcf.gz");
    Map<String, FeatureReader<VariantContext>> variantReaders = new HashMap<>();

    FeatureCodec<VariantContext, ?> codec = new VCFCodec();
    FeatureReader<VariantContext> reader_t6 =
      AbstractFeatureReader.getFeatureReader(t6.getAbsolutePath(), codec, false);
    FeatureReader<VariantContext> reader_t7 =
      AbstractFeatureReader.getFeatureReader(t7.getAbsolutePath(), codec, false);
    FeatureReader<VariantContext> reader_t8 =
      AbstractFeatureReader.getFeatureReader(t8.getAbsolutePath(), codec, false);
    variantReaders.put(((VCFHeader) reader_t6.getHeader()).getGenotypeSamples().get(0), reader_t6);
//    variantReaders.put(((VCFHeader) reader_t7.getHeader()).getGenotypeSamples().get(0), reader_t7);
//    variantReaders.put(((VCFHeader) reader_t8.getHeader()).getGenotypeSamples().get(0), reader_t8);


    GenomicsDBImportConfiguration.ImportConfiguration importConfiguration =
      generateTestLoaderConfiguration();

    ChromosomeInterval chromosomeInterval =
      new ChromosomeInterval("1", 1, 249250619);
    List<VCFHeader> headers = new ArrayList<>();
    for (Map.Entry<String, FeatureReader<VariantContext>> variant : variantReaders.entrySet()) {
      headers.add((VCFHeader) variant.getValue().getHeader());
    }
    Set<VCFHeaderLine> mergedHeader = VCFUtils.smartMergeHeaders(headers, true);
    GenomicsDBImporter importer = new GenomicsDBImporter(
      variantReaders, mergedHeader, chromosomeInterval, importConfiguration);

    importer.importBatch();
    Assert.assertEquals(importer.isDone(), true);
  }

  private GenomicsDBImportConfiguration.ImportConfiguration generateTestLoaderConfiguration() {
    List<GenomicsDBImportConfiguration.Partition> partitions = new ArrayList<>(2);

    GenomicsDBImportConfiguration.Partition.Builder partition0 =
      GenomicsDBImportConfiguration.Partition.newBuilder();
    GenomicsDBImportConfiguration.Partition p0 =
      partition0
        .setBegin(0)
        .setVcfFileName("junk0")
        .setTiledbWorkspace(TILEDB_WORKSPACE)
        .setTiledbArrayName(ARRAY_FOR_PARTITION0)
        .build();
    GenomicsDBImportConfiguration.Partition.Builder partition1 =
      GenomicsDBImportConfiguration.Partition.newBuilder();
    GenomicsDBImportConfiguration.Partition p1 =
      partition1
        .setBegin(1000000)
        .setVcfFileName("junk1")
        .setTiledbWorkspace(TILEDB_WORKSPACE)
        .setTiledbArrayName(ARRAY_FOR_PARTITION1)
        .build();

    partitions.add(p0);
    partitions.add(p1);

    GenomicsDBImportConfiguration.ImportConfiguration.Builder configBuilder =
      GenomicsDBImportConfiguration.ImportConfiguration.newBuilder();
    GenomicsDBImportConfiguration.ImportConfiguration importConfiguration =
      configBuilder
        .setDeleteAndCreateTiledbArray(true)
        .setDoPingPongBuffering(true)
        .setProduceTiledbArray(true)
        .setNumParallelVcfFiles(1)
        .setSizePerColumnPartition(10000)
        .setRowBasedPartitioning(false)
        .addAllColumnPartitions(partitions)
        .build();

    return importConfiguration;
  }
}
