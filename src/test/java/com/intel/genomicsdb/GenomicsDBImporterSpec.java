package com.intel.genomicsdb;

import htsjdk.variant.vcf.VCFFileReader;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class GenomicsDBImporterSpec {

  private static final String TILEDB_WORKSPACE = "./__workspace";
  private static final String ARRAY_FOR_PARTITION0 = "array_test0";
  private static final String ARRAY_FOR_PARTITION1 = "array_test1";

  @Test(groups = {"genomicsdb importer with an interval and multiple GVCFs"})
  public void testMultiGVCFInputs() throws IOException {

    File t6 = new File("tests/inputs/vcfs/t6.vcf.gz");
    File t7 = new File("tests/inputs/vcfs/t7.vcf.gz");
    File t8 = new File("tests/inputs/vcfs/t8.vcf.gz");
    List< VCFFileReader> variantReaders = new ArrayList<>();

    VCFFileReader reader_t6 = new VCFFileReader(t6);
    VCFFileReader reader_t7 = new VCFFileReader(t7);
    VCFFileReader reader_t8 = new VCFFileReader(t8);
    variantReaders.add(reader_t6);
    variantReaders.add(reader_t7);
    variantReaders.add(reader_t8);


    GenomicsDBImportConfiguration.ImportConfiguration importConfiguration =
      generateTestLoaderConfiguration();

    ChromosomeInterval chromosomeInterval =
      new ChromosomeInterval("1",1, 249250621);
    GenomicsDBImporter importer =
      new GenomicsDBImporter(variantReaders, chromosomeInterval, null);

    importer.importBatch();
    assert importer.isDone();
  }

  private GenomicsDBImportConfiguration.ImportConfiguration generateTestLoaderConfiguration() {
    List<GenomicsDBImportConfiguration.Partition> partitions = new ArrayList<>(2);

    GenomicsDBImportConfiguration.TileDBConfig.Builder tB0 =
      GenomicsDBImportConfiguration.TileDBConfig.newBuilder();
    GenomicsDBImportConfiguration.TileDBConfig tileDBConfig_part0 =
      tB0
        .setTiledbWorkspace(TILEDB_WORKSPACE)
        .setTiledbArrayName(ARRAY_FOR_PARTITION0)
        .build();
    GenomicsDBImportConfiguration.TileDBConfig.Builder tB1 =
      GenomicsDBImportConfiguration.TileDBConfig.newBuilder();
    GenomicsDBImportConfiguration.TileDBConfig tileDBConfig_part1 =
      tB1
        .setTiledbWorkspace(TILEDB_WORKSPACE)
        .setTiledbArrayName(ARRAY_FOR_PARTITION1)
        .build();
    GenomicsDBImportConfiguration.Partition.Builder partition0 =
      GenomicsDBImportConfiguration.Partition.newBuilder();
    GenomicsDBImportConfiguration.Partition p0 =
      partition0
        .setBegin(0)
        .setVcfFileName("junk0")
        .setTiledbConfig(tileDBConfig_part0)
        .build();
    GenomicsDBImportConfiguration.Partition.Builder partition1 =
      GenomicsDBImportConfiguration.Partition.newBuilder();
    GenomicsDBImportConfiguration.Partition p1 =
      partition1
        .setBegin(1000000)
        .setVcfFileName("junk1")
        .setTiledbConfig(tileDBConfig_part1)
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
