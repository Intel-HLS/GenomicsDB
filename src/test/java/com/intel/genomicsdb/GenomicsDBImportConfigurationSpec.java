package com.intel.genomicsdb;

import org.testng.annotations.Test;
import com.intel.genomicsdb.GenomicsDBImportConfiguration.ImportConfiguration;
import com.intel.genomicsdb.GenomicsDBImportConfiguration.TileDBConfig;
import com.intel.genomicsdb.GenomicsDBImportConfiguration.Partition;

import java.util.ArrayList;
import java.util.List;

public class GenomicsDBImportConfigurationSpec {

  @Test(groups = {"configuration tests"})
  public void testImportConfiguration() {
    List<Partition.Builder> partitions = new ArrayList<>(2);
    for (Partition.Builder p : partitions) {
      p = Partition.newBuilder();
    }
    ImportConfiguration.Builder importConfiguration = ImportConfiguration.newBuilder();
    importConfiguration.setDeleteAndCreateTiledbArray(true)
      .setDoPingPongBuffering(true)
      .setProduceTiledbArray(true)
      .setNumParallelVcfFiles(1)
      .setSizePerColumnPartition(1000);
    int index = 0;
    for (Partition.Builder p : partitions) {
      importConfiguration.setColumnPartitions(index++,p);
    }
    importConfiguration.build();

    assert importConfiguration.isInitialized();
    assert !importConfiguration.hasCallsetMappingFile();
    assert importConfiguration.hasDoPingPongBuffering();
    assert importConfiguration.getDoPingPongBuffering();

//    System.out.println(importConfiguration.toString());
  }
}
