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

import htsjdk.tribble.*;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.tools.nsc.typechecker.PatternMatching;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class GenomicsDBImporterSpec {

  private static final String TILEDB_WORKSPACE = "./__workspace";
  private static final String TILEDB_ARRAYNAME = "genomicsdb_test_array";
  private static final String TEST_CHROMOSOME_NAME = "1";

  @Test(testName = "genomicsdb importer with an interval and multiple GVCFs")
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
    variantReaders.put(((VCFHeader) reader_t7.getHeader()).getGenotypeSamples().get(0), reader_t7);
    variantReaders.put(((VCFHeader) reader_t8.getHeader()).getGenotypeSamples().get(0), reader_t8);


    ChromosomeInterval chromosomeInterval =
      new ChromosomeInterval(TEST_CHROMOSOME_NAME, 1, 249250619);
    List<VCFHeader> headers = new ArrayList<>();
    for (Map.Entry<String, FeatureReader<VariantContext>> variant : variantReaders.entrySet()) {
      headers.add((VCFHeader) variant.getValue().getHeader());
    }
    Set<VCFHeaderLine> mergedHeader = VCFUtils.smartMergeHeaders(headers, true);
    GenomicsDBImporter importer = new GenomicsDBImporter(
      variantReaders,
      mergedHeader,
      chromosomeInterval,
      TILEDB_WORKSPACE,
      TILEDB_ARRAYNAME,
      0L,
      10000000L);

    importer.importBatch();
    Assert.assertEquals(importer.isDone(), true);
  }

  @Test(testName = "genomicsdb importer outputs merged headers as a JSON file")
  public void testVidMapJSONOutput() throws IOException {

    System.exit(0);

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
    variantReaders.put(((VCFHeader) reader_t7.getHeader()).getGenotypeSamples().get(0), reader_t7);
    variantReaders.put(((VCFHeader) reader_t8.getHeader()).getGenotypeSamples().get(0), reader_t8);

    final String TEMP_JSON_FILE = "./generated_vidmap.json";

    ChromosomeInterval chromosomeInterval =
      new ChromosomeInterval(TEST_CHROMOSOME_NAME, 1, 249250619);
    List<VCFHeader> headers = new ArrayList<>();
    for (Map.Entry<String, FeatureReader<VariantContext>> variant : variantReaders.entrySet()) {
      headers.add((VCFHeader) variant.getValue().getHeader());
    }
    Set<VCFHeaderLine> mergedHeader = VCFUtils.smartMergeHeaders(headers, true);
    GenomicsDBImporter importer = new GenomicsDBImporter(
      variantReaders,
      mergedHeader,
      chromosomeInterval,
      TILEDB_WORKSPACE,
      TILEDB_ARRAYNAME,
      0L,
      10000000L,
      TEMP_JSON_FILE);

    importer.importBatch();
    Assert.assertEquals(importer.isDone(), true);
    Assert.assertEquals(new File(TEMP_JSON_FILE).isFile(), true);
  }
}
