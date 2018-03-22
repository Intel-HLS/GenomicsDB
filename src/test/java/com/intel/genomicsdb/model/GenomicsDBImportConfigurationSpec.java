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
import java.util.List;

public class GenomicsDBImportConfigurationSpec {

  private static final String ARRAY_FOR_PARTITION0 = "array0";
  private static final String ARRAY_FOR_PARTITION1 = "array1";
  private final String TILEDB_WORKSPACE = "/path/to/junk/folder";

  @Test(groups = {"configuration tests"})
  public void testImportConfiguration() {
    List<GenomicsDBImportConfiguration.Partition> partitions = new ArrayList<>(2);

    Coordinates.GenomicsDBColumn genomicsDBColumn0 = Coordinates.GenomicsDBColumn.newBuilder().setTiledbColumn(0).build();
    GenomicsDBImportConfiguration.Partition.Builder partition0 = GenomicsDBImportConfiguration.Partition.newBuilder();
    GenomicsDBImportConfiguration.Partition p0 = partition0.setVcfOutputFilename("junk0").setWorkspace(TILEDB_WORKSPACE)
            .setArray(ARRAY_FOR_PARTITION0).setBegin(genomicsDBColumn0).build();

    Coordinates.GenomicsDBColumn genomicsDBColumn1 = Coordinates.GenomicsDBColumn.newBuilder().setTiledbColumn(0).build();
    GenomicsDBImportConfiguration.Partition.Builder partition1 = GenomicsDBImportConfiguration.Partition.newBuilder();
    GenomicsDBImportConfiguration.Partition p1 = partition1.setVcfOutputFilename("junk1").setWorkspace(TILEDB_WORKSPACE)
            .setArray(ARRAY_FOR_PARTITION1).setBegin(genomicsDBColumn1).build();

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
        .setFailIfUpdating(false)
        .setDisableSyncedWrites(true)
        .setIgnoreCellsNotInPartition(false)
        .build();
    
    Assert.assertEquals(importConfiguration.isInitialized(), true);

    // Assert has methods
    Assert.assertEquals(importConfiguration.hasCallsetMappingFile(), false);
    Assert.assertEquals(importConfiguration.hasDoPingPongBuffering(), true);
    Assert.assertEquals(importConfiguration.hasDeleteAndCreateTiledbArray(), true);
    Assert.assertEquals(importConfiguration.hasDiscardVcfIndex(), false);
    Assert.assertEquals(importConfiguration.hasNumParallelVcfFiles(), true);
    Assert.assertEquals(importConfiguration.hasOffloadVcfOutputProcessing(), false);
    Assert.assertEquals(importConfiguration.hasProduceCombinedVcf(), false);
    Assert.assertEquals(importConfiguration.hasProduceTiledbArray(), true);
    Assert.assertEquals(importConfiguration.hasSizePerColumnPartition(), true);
    Assert.assertEquals(importConfiguration.hasDisableSyncedWrites(), true);
    Assert.assertEquals(importConfiguration.hasIgnoreCellsNotInPartition(), true);

    // Assert gets
    Assert.assertEquals(importConfiguration.getDoPingPongBuffering(), true);
    Assert.assertEquals(importConfiguration.getDoPingPongBuffering(), true);
    Assert.assertEquals(importConfiguration.getDeleteAndCreateTiledbArray(), true);
    Assert.assertEquals(importConfiguration.getNumParallelVcfFiles(), 1);
    Assert.assertEquals(importConfiguration.getProduceTiledbArray(), true);
    Assert.assertEquals(importConfiguration.getSizePerColumnPartition(), 10000);
    Assert.assertEquals(importConfiguration.getColumnPartitionsCount(), 2);
    Assert.assertEquals(importConfiguration.getDisableSyncedWrites(), true);
    Assert.assertEquals(importConfiguration.getIgnoreCellsNotInPartition(), false);
    Assert.assertSame(importConfiguration.getFailIfUpdating(), false);
  }
}
