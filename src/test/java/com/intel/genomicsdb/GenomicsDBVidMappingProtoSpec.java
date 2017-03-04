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

import org.testng.Assert;
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
        .addVcfFieldClassType("FORMAT")
        .build();
    GenomicsDBVidMapProto.InfoField.Builder infoField3 =
      GenomicsDBVidMapProto.InfoField.newBuilder();
    GenomicsDBVidMapProto.InfoField gt =
      infoField3
        .setName("PASS")
        .setType("int")
        .setLength("P")
        .addVcfFieldClassType("FORMAT")
        .build();
    GenomicsDBVidMapProto.InfoField.Builder infoFormatField =
      GenomicsDBVidMapProto.InfoField.newBuilder();
    GenomicsDBVidMapProto.InfoField dp =
      infoFormatField
        .setName("DP")
        .setType("int")
        .setLength("R")
        .addVcfFieldClassType("FORMAT")
        .addVcfFieldClassType("INFO")
        .build();
    infoFields.add(pass);
    infoFields.add(lowQual);
    infoFields.add(pgt);
    infoFields.add(gt);
    infoFields.add(dp);

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

    Assert.assertEquals(vidmap.isInitialized(), true);

    // Assert info fields
    Assert.assertEquals(vidmap.getInfofieldsCount(), 5);
    Assert.assertSame(vidmap.getInfofields(0), pass);
    Assert.assertSame(vidmap.getInfofields(1), lowQual);
    Assert.assertSame(vidmap.getInfofields(2), pgt);
    Assert.assertSame(vidmap.getInfofields(3), gt);
    Assert.assertSame(vidmap.getInfofields(4), dp);

    // Assert chromosomes
    Assert.assertEquals(vidmap.getChromosomesCount(), 4);
    Assert.assertSame(vidmap.getChromosomes(0), one);
    Assert.assertSame(vidmap.getChromosomes(1), two);
    Assert.assertSame(vidmap.getChromosomes(2), three);
    Assert.assertSame(vidmap.getChromosomes(3), gl000195);
  }
}



