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

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.testng.annotations.DataProvider;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public final class GenomicsDBTestUtils {

    @DataProvider(name = "vcfFile")
    public static Object[][] vcfFile() {
        File t0 = new File("tests/inputs/vcfs/t0.vcf.gz");

        Map<String, FeatureReader<VariantContext>> sampleToReaderMap = new HashMap<>();
        FeatureReader<VariantContext> reader_t0 =
                AbstractFeatureReader.getFeatureReader(t0.getAbsolutePath(), new VCFCodec(), false);
        sampleToReaderMap.put(((VCFHeader) reader_t0.getHeader()).getGenotypeSamples().get(0), reader_t0);

        return new Object[][]{{sampleToReaderMap}};
    }

    @DataProvider(name = "vcfFiles")
    public static Object[][] vcfFiles() {
        File t6 = new File("tests/inputs/vcfs/t6.vcf.gz");
        File t7 = new File("tests/inputs/vcfs/t7.vcf.gz");
        File t8 = new File("tests/inputs/vcfs/t8.vcf.gz");

        Map<String, FeatureReader<VariantContext>> sampleToReaderMap = new HashMap<>();
        FeatureReader<VariantContext> reader_t6 =
                AbstractFeatureReader.getFeatureReader(t6.getAbsolutePath(), new VCFCodec(), false);
        FeatureReader<VariantContext> reader_t7 =
                AbstractFeatureReader.getFeatureReader(t7.getAbsolutePath(), new VCFCodec(), false);
        FeatureReader<VariantContext> reader_t8 =
                AbstractFeatureReader.getFeatureReader(t8.getAbsolutePath(), new VCFCodec(), false);
        sampleToReaderMap.put(((VCFHeader) reader_t6.getHeader()).getGenotypeSamples().get(0), reader_t6);
        sampleToReaderMap.put(((VCFHeader) reader_t7.getHeader()).getGenotypeSamples().get(0), reader_t7);
        sampleToReaderMap.put(((VCFHeader) reader_t8.getHeader()).getGenotypeSamples().get(0), reader_t8);

        return new Object[][]{{sampleToReaderMap}};
    }

    @DataProvider(name = "nullFeatureReaders")
    public static Object[][] nullFeatureReaders() {
        Map<String, FeatureReader<VariantContext>> sampleToReaderMap = new HashMap<>();
        sampleToReaderMap.put("ABC", null);

        return new Object[][]{{sampleToReaderMap}};
    }

  @DataProvider(name="loaderQueryHostFilesTest0")
  public static Object[][] loaderQueryHostFilesTest0() throws IOException {
    String queryPreamble = "{\n\"workspace\": \"hdfs://tmp/ws\",\n\"array\": \"part\",\n";
    String loaderPreamble = "{\n\"row_based_partitioning\": false,\n\"column_partitions\" : [\n";
    String loaderEpilogue = "\n],\n\"vid_mapping_file\" : \"tests/inputs/vid.json\"\n}";
    String loaderPartEpilogue = "\"workspace\": \"hdfs://tmp/ws\", \"vcf_output_filename\": \"/tmp/test0.vcf.gz\"}";

    String queryCol = "\"query_column_ranges\" : [ [ ";
    String loaderPart = "";
    for(int i=0; i<3; i++) {
       queryCol += "[ "+ String.valueOf(i*500000+1000) + "," + String.valueOf(i*500000+1100) + "]";
       loaderPart += "{\"begin\": "+String.valueOf(i*500000)+", \"array\": \"part"+i+"\","+loaderPartEpilogue;
       if (i<2) {
         queryCol += ",";
         loaderPart += ",\n";
       }
    }
    String query = queryPreamble + queryCol + " ] ]\n}"; 
    String loader = loaderPreamble + loaderPart + loaderEpilogue;
    File tmpQFile = File.createTempFile("query", ".json");
    tmpQFile.deleteOnExit();
    FileWriter fQ = new FileWriter(tmpQFile);
    fQ.write(query);
    fQ.close();
    File tmpLFile = File.createTempFile("loader", ".json");
    tmpLFile.deleteOnExit();
    FileWriter fL = new FileWriter(tmpLFile);
    fL.write(loader);
    fL.close();
    File tmpHFile = File.createTempFile("hostfile", "");
    tmpHFile.deleteOnExit();
    FileWriter fH = new FileWriter(tmpHFile);
    fH.write("localhost");
    fH.close();

    return new Object[][] { { tmpQFile.getAbsolutePath(), tmpLFile.getAbsolutePath(), tmpHFile.getAbsolutePath() } };
  }

  @DataProvider(name="loaderQueryHostFilesTest1")
  public static Object[][] loaderQueryHostFilesTest1() throws IOException {
    String queryPreamble = "{\n\"workspace\": \"hdfs://tmp/ws\",\n\"array\": \"part\",\n";
    String loaderPreamble = "{\n\"row_based_partitioning\": false,\n\"column_partitions\" : [\n";
    String loaderEpilogue = "\n],\n\"vid_mapping_file\" : \"tests/inputs/vid.json\"\n}";
    String loaderPartEpilogue = "\"workspace\": \"hdfs://tmp/ws\", \"vcf_output_filename\": \"/tmp/test0.vcf.gz\"}";

    String queryCol = "\"query_column_ranges\" : [ [ ";
    String loaderPart = "{\"begin\": 0, \"array\": \"part\","+loaderPartEpilogue;
    for(int i=0; i<3; i++) {
       queryCol += "["+ String.valueOf(i*500000+2000) + "," + String.valueOf(i*500000+2100) + "]";
       if (i<2) {
         queryCol += ",";
       }
    }
    String query = queryPreamble + queryCol + " ] ]\n}"; 
    String loader = loaderPreamble + loaderPart + loaderEpilogue;
    File tmpQFile = File.createTempFile("query", ".json");
    tmpQFile.deleteOnExit();
    FileWriter fQ = new FileWriter(tmpQFile);
    fQ.write(query);
    fQ.close();
    File tmpLFile = File.createTempFile("loader", ".json");
    tmpLFile.deleteOnExit();
    FileWriter fL = new FileWriter(tmpLFile);
    fL.write(loader);
    fL.close();
    File tmpHFile = File.createTempFile("hostfile", "");
    tmpHFile.deleteOnExit();
    FileWriter fH = new FileWriter(tmpHFile);
    fH.write("localhost");
    fH.close();

    return new Object[][] { { tmpQFile.getAbsolutePath(), tmpLFile.getAbsolutePath(), tmpHFile.getAbsolutePath() } };
  }

  @DataProvider(name="loaderQueryHostFilesTest2")
  public static Object[][] loaderQueryHostFilesTest2() throws IOException {
    String queryPreamble = "{\n\"workspace\": \"hdfs://tmp/ws\",\n\"array\": \"part\",\n";
    String loaderPreamble = "{\n\"row_based_partitioning\": false,\n\"column_partitions\" : [\n";
    String loaderEpilogue = "\n],\n\"vid_mapping_file\" : \"tests/inputs/vid.json\"\n}";
    String loaderPartEpilogue = "\"workspace\": \"hdfs://tmp/ws\", \"vcf_output_filename\": \"/tmp/test0.vcf.gz\"}";

    String queryCol = "\"query_column_ranges\" : [ [ [ 500, 25000]";
    String loaderPart = "{\"begin\": 0, \"array\": \"part\","+loaderPartEpilogue;
    String query = queryPreamble + queryCol + " ] ]\n}"; 
    String loader = loaderPreamble + loaderPart + loaderEpilogue;
    File tmpQFile = File.createTempFile("query", ".json");
    tmpQFile.deleteOnExit();
    FileWriter fQ = new FileWriter(tmpQFile);
    fQ.write(query);
    fQ.close();
    File tmpLFile = File.createTempFile("loader", ".json");
    tmpLFile.deleteOnExit();
    FileWriter fL = new FileWriter(tmpLFile);
    fL.write(loader);
    fL.close();
    File tmpHFile = File.createTempFile("hostfile", "");
    tmpHFile.deleteOnExit();
    FileWriter fH = new FileWriter(tmpHFile);
    fH.write("localhost");
    fH.close();

    return new Object[][] { { tmpQFile.getAbsolutePath(), tmpLFile.getAbsolutePath(), tmpHFile.getAbsolutePath() } };
  }

  @DataProvider(name="loaderQueryHostFilesTest3")
  public static Object[][] loaderQueryHostFilesTest3() throws IOException {
    String queryPreamble = "{\n\"workspace\": \"/tmp/ws\",\n\"array\": \"part\",\n";
    String loaderPreamble = "{\n\"row_based_partitioning\": false,\n\"column_partitions\" : [\n";
    String loaderEpilogue = "\n],\n\"vid_mapping_file\" : \"tests/inputs/vid.json\"\n}";
    String loaderPartEpilogue = "\"workspace\": \"/tmp/ws\", \"vcf_output_filename\": \"/tmp/test0.vcf.gz\"}";

    String queryCol = "\"query_column_ranges\" : [ [ [ 500, 25000]";
    String loaderPart = "{\"begin\": 0, \"array\": \"part\","+loaderPartEpilogue;
    String query = queryPreamble + queryCol + " ] ]\n}"; 
    String loader = loaderPreamble + loaderPart + loaderEpilogue;
    File tmpQFile = File.createTempFile("query", ".json");
    tmpQFile.deleteOnExit();
    FileWriter fQ = new FileWriter(tmpQFile);
    fQ.write(query);
    fQ.close();
    File tmpLFile = File.createTempFile("loader", ".json");
    tmpLFile.deleteOnExit();
    FileWriter fL = new FileWriter(tmpLFile);
    fL.write(loader);
    fL.close();
    File tmpHFile = File.createTempFile("hostfile", "");
    tmpHFile.deleteOnExit();
    FileWriter fH = new FileWriter(tmpHFile);
    fH.write("node0\nnode1\nnode2\nnode3\n");
    fH.close();

    return new Object[][] { { tmpQFile.getAbsolutePath(), tmpLFile.getAbsolutePath(), tmpHFile.getAbsolutePath() } };
  }
}
