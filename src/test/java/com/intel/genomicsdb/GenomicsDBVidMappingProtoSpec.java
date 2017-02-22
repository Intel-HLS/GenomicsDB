package com.intel.genomicsdb;

import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class GenomicsDBVidMappingProtoSpec {

  @Test(groups = {"vid mapping tests"})
  public void testGenomicsDBVidMapProto() {

    // Create INFO fields
    List<GenomicsDBVidMapProto.InfoField> infoFields = new ArrayList<>(4);
    GenomicsDBVidMapProto.InfoField.Builder infoField0 =
      GenomicsDBVidMapProto.InfoField.newBuilder();
    GenomicsDBVidMapProto.InfoField pass =
      infoField0
        .setName("PASS")
        .setType("int")
        .build();
    GenomicsDBVidMapProto.InfoField.Builder infoField1 =
      GenomicsDBVidMapProto.InfoField.newBuilder();
    GenomicsDBVidMapProto.InfoField lowQual =
      infoField1
        .setName("LowQual")
        .setType("int")
        .build();
    GenomicsDBVidMapProto.InfoField.Builder infoField2 =
      GenomicsDBVidMapProto.InfoField.newBuilder();
    GenomicsDBVidMapProto.InfoField pgt =
      infoField2
        .setName("PGT")
        .setType("char")
        .setLength("VAR")
        .setVcfFieldClassType("FORMAT")
        .build();
    GenomicsDBVidMapProto.InfoField.Builder infoField3 =
      GenomicsDBVidMapProto.InfoField.newBuilder();
    GenomicsDBVidMapProto.InfoField gt =
      infoField3
        .setName("PASS")
        .setType("int")
        .setLength("P")
        .setVcfFieldClassType("FORMAT")
        .build();
    infoFields.add(pass);
    infoFields.add(lowQual);
    infoFields.add(pgt);
    infoFields.add(gt);

    // Create chromosomes
    List<GenomicsDBVidMapProto.Chromosome> chromosomes = new ArrayList<>(4);
    GenomicsDBVidMapProto.Chromosome.Builder contig0 =
      GenomicsDBVidMapProto.Chromosome.newBuilder();
    GenomicsDBVidMapProto.Chromosome one =
      contig0.setName("1").setLength(249250621).setTiledbColumnOffset(0).build();
    GenomicsDBVidMapProto.Chromosome.Builder contig1 =
      GenomicsDBVidMapProto.Chromosome.newBuilder();
    GenomicsDBVidMapProto.Chromosome two =
      contig1
        .setName("2")
        .setLength(243199373)
        .setTiledbColumnOffset(249250621)
        .build();
    GenomicsDBVidMapProto.Chromosome.Builder contig2 =
      GenomicsDBVidMapProto.Chromosome.newBuilder();
    GenomicsDBVidMapProto.Chromosome three =
      contig2
        .setName("3")
        .setLength(243199373)
        .setTiledbColumnOffset(249250621)
        .build();
    GenomicsDBVidMapProto.Chromosome.Builder contig3 =
      GenomicsDBVidMapProto.Chromosome.newBuilder();
    GenomicsDBVidMapProto.Chromosome gl000195 =
      contig3
        .setName("GL000195.1")
        .setLength(182896)
        .setTiledbColumnOffset(3099921162L)
        .build();
    chromosomes.add(one);
    chromosomes.add(two);
    chromosomes.add(three);
    chromosomes.add(gl000195);

    GenomicsDBVidMapProto.VidMapping.Builder vidmapBuilder =
      GenomicsDBVidMapProto.VidMapping.newBuilder();
    GenomicsDBVidMapProto.VidMapping vidmap =
      vidmapBuilder
        .addAllChromosomes(chromosomes)
        .addAllInfofields(infoFields)
        .build();

    assert vidmap.isInitialized();

    // Assert info fields
    assert vidmap.getInfofieldsCount() == 4;
    assert vidmap.getInfofields(0) == pass;
    assert vidmap.getInfofields(1) == lowQual;
    assert vidmap.getInfofields(2) == pgt;
    assert vidmap.getInfofields(3) == gt;

    // Assert chromosomes
    assert vidmap.getChromosomesCount() == 4;
    assert vidmap.getChromosomes(0) == one;
    assert vidmap.getChromosomes(1) == two;
    assert vidmap.getChromosomes(2) == three;
    assert vidmap.getChromosomes(3) == gl000195;
  }
}



