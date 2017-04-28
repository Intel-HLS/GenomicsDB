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

import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.apache.commons.io.FileUtils;
import org.junit.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static com.googlecode.protobuf.format.JsonFormat.printToString;

public final class GenomicsDBFeatureReaderSpec {

  private static File TEST_TILEDB_WORKSPACE = new File("__workspace");
  private static String TEST_TILEDB_ARRAY_NAME = "featureReaderTest";
  private static File TEMP_VIDMAP_JSON_FILE = new File("generated_vid_map.json");
  private static File TEMP_CALLSETMAP_JSON_FILE = new File("generated_callset_map.json");
  private static File TEST_LOADER_JSON_FILE = new File("generated_loader.json");
  private static File TEST_QUERY_JSON_FILE = new File("generated_query.json");
  private static File TEST_REFERENCE_GENOME = new File("tests/inputs/Homo_sapiens_assembly19.fasta");
  private static File TEMP_VID_JSON_FILE = new File("generated_vidmap.json");
  private static File TEMP_CALLSET_JSON_FILE = new File("generated_callsetmap.json");

  private static final String TEST_CHROMOSOME_NAME = "1";

  static File printQueryJSONFile(
    GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration,
    String filename) {
    String queryJSONString = printToString(exportConfiguration);

    File tempQueryJSONFile = new File(filename);

    try( PrintWriter out = new PrintWriter(tempQueryJSONFile)  ){
      out.println(queryJSONString);
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }
    return tempQueryJSONFile;
  }

  @Test(testName = "Feature Reader with Merged Header",
        dataProvider = "vcfFiles",
        dataProviderClass = GenomicsDBTestUtils.class)
  public void testFeatureReaderWithMergedHeader(
    Map<String, FeatureReader<VariantContext>> variantReaders) throws IOException {



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
      TEST_TILEDB_WORKSPACE.getAbsolutePath(),
      TEST_TILEDB_ARRAY_NAME,
      1000L,
      10000000L,
      TEMP_VID_JSON_FILE.getAbsolutePath(),
      TEMP_CALLSET_JSON_FILE.getAbsolutePath());

    importer.importBatch();

    GenomicsDBImportConfiguration.ImportConfiguration.Builder builder =
      GenomicsDBImportConfiguration.ImportConfiguration.newBuilder();

    GenomicsDBImportConfiguration.Partition.Builder pBuilder =
      GenomicsDBImportConfiguration.Partition.newBuilder();

    GenomicsDBImportConfiguration.Partition partition =
      pBuilder
      .setArray(TEST_TILEDB_ARRAY_NAME)
      .setWorkspace(TEST_TILEDB_WORKSPACE.getAbsolutePath())
      .setBegin(0)
      .build();

    List<GenomicsDBImportConfiguration.Partition> partitionList = new ArrayList<>();
    partitionList.add(partition);

    GenomicsDBImportConfiguration.ImportConfiguration importConfiguration =
      builder
        .setProduceTiledbArray(true)
        .setCallsetMappingFile(TEMP_CALLSETMAP_JSON_FILE.getAbsolutePath())
        .setVidMappingFile(TEMP_VIDMAP_JSON_FILE.getAbsolutePath())
        .setSizePerColumnPartition(1000L)
        .setSegmentSize(1048576)
        .addAllColumnPartitions(partitionList)
        .build();

    File importJSONFile = GenomicsDBImporter.printLoaderJSONFile(
      importConfiguration, TEST_LOADER_JSON_FILE.getAbsolutePath());

    GenomicsDBExportConfiguration.ExportConfiguration.Builder eBuilder =
      GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();

    GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration =
      eBuilder
        .setWorkspace(TEST_TILEDB_WORKSPACE.getAbsolutePath())
        .setArray(TEST_TILEDB_ARRAY_NAME)
        .setReferenceGenome(TEST_REFERENCE_GENOME.getAbsolutePath())
        .build();

//    File queryJSONFile = printQueryJSONFile(exportConfiguration, TEST_QUERY_JSON_FILE.getAbsolutePath());
//    BCF2Codec codec = new BCF2Codec();
//    GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> reader =
//      new GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream>(
//        importJSONFile.getAbsolutePath(), queryJSONFile.getAbsolutePath(), codec);

    FileUtils.deleteQuietly(importJSONFile);
  }
}
