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
    GenomicsDBVidMapProto.FieldLengthDescriptorComponentPB.Builder lengthDescriptorComponentBuilder =
        GenomicsDBVidMapProto.FieldLengthDescriptorComponentPB.newBuilder();
    lengthDescriptorComponentBuilder.setVariableLengthDescriptor("VAR");
    GenomicsDBVidMapProto.InfoField pgt =
      infoField2
        .setName("PGT")
        .setType("char")
        .addLength(lengthDescriptorComponentBuilder.build())
        .addVcfFieldClass("FORMAT")
        .build();
    GenomicsDBVidMapProto.InfoField.Builder infoField3 =
      GenomicsDBVidMapProto.InfoField.newBuilder();
    lengthDescriptorComponentBuilder.setVariableLengthDescriptor("P");
    GenomicsDBVidMapProto.InfoField gt =
      infoField3
        .setName("PASS")
        .setType("int")
        .addLength(lengthDescriptorComponentBuilder.build())
        .addVcfFieldClass("FORMAT")
        .build();
    GenomicsDBVidMapProto.InfoField.Builder infoFormatField =
      GenomicsDBVidMapProto.InfoField.newBuilder();
    lengthDescriptorComponentBuilder.setVariableLengthDescriptor("R");
    GenomicsDBVidMapProto.InfoField dp =
      infoFormatField
        .setName("DP")
        .setType("int")
        .addLength(lengthDescriptorComponentBuilder.build())
        .addVcfFieldClass("FORMAT")
        .addVcfFieldClass("INFO")
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

    GenomicsDBVidMapProto.VidMappingPB.Builder vidmapBuilder =
      GenomicsDBVidMapProto.VidMappingPB.newBuilder();
    GenomicsDBVidMapProto.VidMappingPB vidmap =
      vidmapBuilder
        .addAllContigs(chromosomes)
        .addAllFields(infoFields)
        .build();

    Assert.assertEquals(vidmap.isInitialized(), true);

    // Assert info fields
    Assert.assertEquals(vidmap.getFieldsCount(), 5);
    Assert.assertSame(vidmap.getFields(0), pass);
    Assert.assertSame(vidmap.getFields(1), lowQual);
    Assert.assertSame(vidmap.getFields(2), pgt);
    Assert.assertSame(vidmap.getFields(3), gt);
    Assert.assertSame(vidmap.getFields(4), dp);

    // Assert chromosomes
    Assert.assertEquals(vidmap.getContigsCount(), 4);
    Assert.assertSame(vidmap.getContigs(0), one);
    Assert.assertSame(vidmap.getContigs(1), two);
    Assert.assertSame(vidmap.getContigs(2), three);
    Assert.assertSame(vidmap.getContigs(3), gl000195);
  }
}



