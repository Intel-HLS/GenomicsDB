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

import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public final class GenomicsDBFeatureReaderSpec {

  private static String TEST_TILEDB_WORKSPACE = "./__workspace";
  private static String TEST_TILEDB_ARRAY_NAME = "featureReaderTest";
  private static String TEMP_VIDMAP_JSON_FILE = "./generated_vid_map.json";
  private static String TEMP_CALLSETMAP_JSON_FILE = "./generated_callset_map.json";

  @Test(testName = "Feature Reader with Merged Header")
  public void testFeatureReaderWithMergedHeader() {

    GenomicsDBImportConfiguration.ImportConfiguration.Builder builder =
      GenomicsDBImportConfiguration.ImportConfiguration.newBuilder();

    GenomicsDBImportConfiguration.Partition.Builder pBuilder =
      GenomicsDBImportConfiguration.Partition.newBuilder();

    GenomicsDBImportConfiguration.Partition partition =
      pBuilder
      .setArray(TEST_TILEDB_ARRAY_NAME)
      .setWorkspace(TEST_TILEDB_WORKSPACE)
      .setBegin(0)
      .build();

    List<GenomicsDBImportConfiguration.Partition> partitionList = new ArrayList<>();
    partitionList.add(partition);

    GenomicsDBImportConfiguration.ImportConfiguration importConfiguration =
      builder
        .setProduceTiledbArray(true)
        .setCallsetMappingFile(TEMP_CALLSETMAP_JSON_FILE)
        .setVidMappingFile(TEMP_VIDMAP_JSON_FILE)
        .setSizePerColumnPartition(1000L)
        .setSegmentSize(1048576)
        .addAllColumnPartitions(partitionList)
        .build();

//    GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> reader =
//      GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream>();
  }
}
