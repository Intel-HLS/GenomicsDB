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
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.apache.commons.io.FileUtils;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.*;

public final class GenomicsDBImporterSpec {

  private static final String WORKSPACE = "./__workspace";
  private static final String TILEDB_ARRAYNAME = "genomicsdb_test_array";
  private static final String TEST_CHROMOSOME_NAME = "1";

  @DataProvider(name="vcfFiles")
  public Object[][] vcfFiles() {
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

    return new Object[][] {{ variantReaders }};
  }

  @Test(testName = "genomicsdb importer with an interval and multiple GVCFs",
        dataProvider = "vcfFiles")
  public void testMultiGVCFInputs(Map<String, FeatureReader<VariantContext>> variantReaders)
    throws IOException {

    ChromosomeInterval chromosomeInterval =
      new ChromosomeInterval(TEST_CHROMOSOME_NAME, 1, 249250619);

    Set<VCFHeaderLine> mergedHeader = createMergedHeader(variantReaders);

    GenomicsDBImporter importer = new GenomicsDBImporter(
      variantReaders,
      mergedHeader,
      chromosomeInterval,
      WORKSPACE,
      TILEDB_ARRAYNAME,
      0L,
      10000000L);

    importer.importBatch();
    Assert.assertEquals(importer.isDone(), true);
  }

  private Set<VCFHeaderLine> createMergedHeader(
    Map<String, FeatureReader<VariantContext>> variantReaders) {

    List<VCFHeader> headers = new ArrayList<>();
    for (Map.Entry<String, FeatureReader<VariantContext>> variant : variantReaders.entrySet()) {
      headers.add((VCFHeader) variant.getValue().getHeader());
    }

    return VCFUtils.smartMergeHeaders(headers, true);
  }

  @Test(testName = "genomicsdb importer outputs merged headers as a JSON file",
        dataProvider = "vcfFiles")
  public void testVidMapJSONOutput(Map<String, FeatureReader<VariantContext>> variantReaders)
    throws IOException {

    final String TEMP_VID_JSON_FILE = "./generated_vidmap.json";
    final String TEMP_CALLSET_JSON_FILE = "./generated_callsetmap.json";

    ChromosomeInterval chromosomeInterval =
      new ChromosomeInterval(TEST_CHROMOSOME_NAME, 1, 249250619);

    Set<VCFHeaderLine> mergedHeader = createMergedHeader(variantReaders);

    GenomicsDBImporter importer = new GenomicsDBImporter(
      variantReaders,
      mergedHeader,
      chromosomeInterval,
      WORKSPACE,
      TILEDB_ARRAYNAME,
      0L,
      10000000L,
      TEMP_VID_JSON_FILE,
      TEMP_CALLSET_JSON_FILE);

    GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB_A = 
      GenomicsDBImporter.generateSortedCallSetMap(variantReaders, false);

    importer.importBatch();
    Assert.assertEquals(importer.isDone(), true);
    Assert.assertEquals(new File(TEMP_VID_JSON_FILE).isFile(), true);
    Assert.assertEquals(new File(TEMP_CALLSET_JSON_FILE).isFile(), true);

    JSONParser parser = new JSONParser();

    try {
      FileReader fileReader = new FileReader(TEMP_CALLSET_JSON_FILE);
      JSONObject jsonObject =
        (JSONObject) parser.parse(fileReader);

      JSONArray callsetArray = (JSONArray) jsonObject.get("callsets");

      int index = 0;
      for (Object cObject : callsetArray) {
        GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap sampleIDToTileDBIDMap =
          (GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap) cObject;
        String sampleName = sampleIDToTileDBIDMap.getSampleName();
        Long tiledbRowIndex_B = sampleIDToTileDBIDMap.getRowIdx();
        String stream_name_B = sampleIDToTileDBIDMap.getStreamName();

        Long tileDBRowIndex_A =
          callsetMappingPB_A.getCallsets(index).getRowIdx();
        String stream_name_A =
          callsetMappingPB_A.getCallsets(index).getStreamName();

        Assert.assertEquals(tileDBRowIndex_A, tiledbRowIndex_B);
        Assert.assertEquals(stream_name_A, stream_name_B);
        index++;
      }
    } catch (ParseException p) {
      p.printStackTrace();
    }

    FileUtils.deleteQuietly(new File(TEMP_VID_JSON_FILE));
    FileUtils.deleteQuietly(new File(TEMP_CALLSET_JSON_FILE));
  }

  @AfterTest
  public void deleteWorkspace() throws IOException {
    FileUtils.deleteDirectory(new File(WORKSPACE));
  }
}
