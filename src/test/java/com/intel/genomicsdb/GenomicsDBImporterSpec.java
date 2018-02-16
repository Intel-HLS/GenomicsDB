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

import com.intel.genomicsdb.model.BaseImportConfig;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.VariantContext;
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
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public final class GenomicsDBImporterSpec {

    private static final File WORKSPACE = new File("__workspace");
    private static final String TEST_CHROMOSOME_NAME = "1";
    private static final File TEMP_VID_JSON_FILE = new File("generated_vidmap.json");
    private static final File TEMP_CALLSET_JSON_FILE = new File("generated_callsetmap.json");
    private static final File TEMP_VCF_HEADER_FILE = new File("vcfheader.vcf");
    private static final File TEMP_QUERY_JSON_FILE = new File("query.json");

    @Test(testName = "genomicsdb importer with an interval and multiple GVCFs",
            dataProvider = "vcfFiles",
            dataProviderClass = GenomicsDBTestUtils.class)
    public void testMultiGVCFInputs(Map<String, FeatureReader<VariantContext>> sampleToReaderMap)
            throws IOException {

        ChromosomeInterval chromosomeInterval =
                new ChromosomeInterval(TEST_CHROMOSOME_NAME, 1, 249250619);

        Set<VCFHeaderLine> mergedHeader = createMergedHeader(sampleToReaderMap);

        GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB =
                GenomicsDBImporter.generateSortedCallSetMap(sampleToReaderMap, true, false);

        GenomicsDBImporter importer = new GenomicsDBImporter(
                sampleToReaderMap,
                mergedHeader,
                chromosomeInterval,
                WORKSPACE.getAbsolutePath(),
                TEST_CHROMOSOME_NAME,
                1024L,
                10000000L);

        importer.importBatch();

        GenomicsDBImporter.writeCallsetMapJSONFile(TEMP_CALLSET_JSON_FILE.getAbsolutePath(),
                callsetMappingPB);
        Assert.assertEquals(importer.isDone(), true);
        File callsetFile = new File(TEMP_CALLSET_JSON_FILE.getAbsolutePath());
        Assert.assertTrue(callsetFile.exists());
    }

    @Test(testName = "genomicsdb importer outputs merged headers as a JSON file",
            dataProvider = "vcfFiles",
            dataProviderClass = GenomicsDBTestUtils.class)
    public void testVidMapJSONOutput(Map<String, FeatureReader<VariantContext>> sampleToReaderMap)
            throws IOException {

        ChromosomeInterval chromosomeInterval =
                new ChromosomeInterval(TEST_CHROMOSOME_NAME, 1, 249250619);

        Set<VCFHeaderLine> mergedHeaderLines = createMergedHeader(sampleToReaderMap);

        GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB_A =
                GenomicsDBImporter.generateSortedCallSetMap(sampleToReaderMap, true, false);

        GenomicsDBImporter importer = new GenomicsDBImporter(
                sampleToReaderMap,
                mergedHeaderLines,
                chromosomeInterval,
                WORKSPACE.getAbsolutePath(),
                TEST_CHROMOSOME_NAME,
                1024L,
                10000000L);

        GenomicsDBImporter.writeVidMapJSONFile(TEMP_VID_JSON_FILE.getAbsolutePath(), mergedHeaderLines);
        GenomicsDBImporter.writeVcfHeaderFile(TEMP_VCF_HEADER_FILE.getAbsolutePath(), mergedHeaderLines);
        GenomicsDBImporter.writeCallsetMapJSONFile(TEMP_CALLSET_JSON_FILE.getAbsolutePath(),
                callsetMappingPB_A);

        importer.importBatch();
        Assert.assertEquals(importer.isDone(), true);
        Assert.assertEquals(TEMP_VID_JSON_FILE.isFile(), true);
        Assert.assertEquals(TEMP_CALLSET_JSON_FILE.isFile(), true);
        Assert.assertEquals(TEMP_VCF_HEADER_FILE.isFile(), true);

        JSONParser parser = new JSONParser();

        try {
            FileReader fileReader = new FileReader(TEMP_CALLSET_JSON_FILE.getAbsolutePath());
            JSONObject jsonObject = (JSONObject) parser.parse(fileReader);
            JSONArray callsetArray = (JSONArray) jsonObject.get("callsets");

            int index = 0;

            for (Object cObject : callsetArray) {
                JSONObject sampleObject = (JSONObject) cObject;
                String sampleName_B = (String) sampleObject.get("sample_name");
                Long tiledbRowIndex_B = (Long) sampleObject.get("row_idx");
                String stream_name_B = (String) sampleObject.get("stream_name");

                String sampleName_A = callsetMappingPB_A.getCallsets(index).getSampleName();
                Long tileDBRowIndex_A = callsetMappingPB_A.getCallsets(index).getRowIdx();
                String stream_name_A = callsetMappingPB_A.getCallsets(index).getStreamName();

                Assert.assertEquals(sampleName_A, sampleName_B);
                Assert.assertEquals(tileDBRowIndex_A, tiledbRowIndex_B);
                Assert.assertEquals(stream_name_A, stream_name_B);
                index++;
            }
        } catch (ParseException p) {
            p.printStackTrace();
        }
    }

    @Test(testName = "genomicsdb importer with null feature readers",
            dataProvider = "nullFeatureReaders",
            dataProviderClass = GenomicsDBTestUtils.class,
            expectedExceptions = IllegalArgumentException.class)
    public void testNullFeatureReaders(Map<String, FeatureReader<VariantContext>> sampleToReaderMap)
            throws IOException {

        new ChromosomeInterval(TEST_CHROMOSOME_NAME, 1, 249250619);

        GenomicsDBImporter.generateSortedCallSetMap(sampleToReaderMap, true, false);
    }

    @Test(testName = "genomicsdb command line config and parallel import with multiple chromosome intervals and one GVCFs")
    public void testImportUsingCommandLineConfigAndParallelImport() throws IOException {
        long start = System.nanoTime();

        //Given
        String[] args = ("-L 1:12000-13000 -L 1:17000-18000 " +
                "-w " + WORKSPACE.getAbsolutePath() + " " +
                "--vidmap-output " + TEMP_VID_JSON_FILE.getAbsolutePath() + " " +
                "--callset-output " + TEMP_CALLSET_JSON_FILE.getAbsolutePath() + " " +
                "tests/inputs/vcfs/t0.vcf.gz").split(" ");
        BaseImportConfig config = new BaseImportConfig("TestGenomicsDBImporterWithMergedVCFHeader", args);
        config.setVcfBufferSizePerColumnPartition(10000L);
        config.setSegmentSize(1048576L);

        //When
        GenomicsDBImporter.parallelImport(config);

        long duration = (System.nanoTime() - start) / 1_000_000;
        System.out.printf("Processed %d chromosomes intervals in %d millis\n", config.getChromosomeIntervalList().size(), duration);

        //Then
        Assert.assertEquals(TEMP_VID_JSON_FILE.isFile(), true);
        Assert.assertEquals(TEMP_CALLSET_JSON_FILE.isFile(), true);

        createQueryJsonFile();
        assert(TEMP_QUERY_JSON_FILE.exists());
        GenomicsDBFeatureReader reader = new GenomicsDBFeatureReader<>("",
                TEMP_QUERY_JSON_FILE.getAbsolutePath(), new BCF2Codec());
        CloseableTribbleIterator<VariantContext> gdbIterator = reader.iterator();
        List<VariantContext> varCtxList = new ArrayList<>();
        while (gdbIterator.hasNext()) varCtxList.add(gdbIterator.next());

        assert(varCtxList.size() == 2);
        assert(varCtxList.get(0).getStart() == 12141);
        assert(varCtxList.get(1).getStart() == 17385);
    }

    @AfterTest
    public void deleteWorkspace() throws IOException {
        if (TEMP_VID_JSON_FILE.exists()) FileUtils.deleteQuietly(TEMP_VID_JSON_FILE);
        if (TEMP_CALLSET_JSON_FILE.exists()) FileUtils.deleteQuietly(TEMP_CALLSET_JSON_FILE);
        if (TEMP_QUERY_JSON_FILE.exists()) FileUtils.deleteQuietly(TEMP_QUERY_JSON_FILE);
        FileUtils.deleteDirectory(WORKSPACE);
    }

    private Set<VCFHeaderLine> createMergedHeader(
            Map<String, FeatureReader<VariantContext>> variantReaders) {

        List<VCFHeader> headers = new ArrayList<>();
        for (Map.Entry<String, FeatureReader<VariantContext>> variant : variantReaders.entrySet()) {
            headers.add((VCFHeader) variant.getValue().getHeader());
        }

        return VCFUtils.smartMergeHeaders(headers, true);
    }

    private String buildQueryJson(String workspace, String arrayName, String callsetMappingFile, String vidMappingFile) {
        JSONObject obj = new JSONObject();
        obj.put("workspace", workspace);
        obj.put("array", arrayName);
        obj.put("scan_full", Boolean.TRUE);
        obj.put("callset_mapping_file", callsetMappingFile);
        obj.put("vid_mapping_file", vidMappingFile);
        JSONArray colrg = new JSONArray();
        JSONArray colrgInn = new JSONArray();
        JSONArray colrgInnInn = new JSONArray();
        colrgInnInn.add(0);
        colrgInnInn.add(10000);
        colrgInn.add(colrgInnInn);
        colrg.add(colrgInn);
        obj.put("xxquery_column_ranges", colrg);
        JSONArray rowrg = new JSONArray();
        JSONArray rowrgInn = new JSONArray();
        rowrgInn.add(0);
        rowrgInn.add(1);
        rowrgInn.add(2);
        rowrg.add(rowrgInn);
        obj.put("xxquery_row_ranges", rowrg);
        obj.put("reference_genome", "tests/inputs/chr1_10MB.fasta.gz");
        obj.put("index_output_VCF", Boolean.TRUE);
        obj.put("produce_GT_field", Boolean.TRUE);
        return obj.toJSONString();
    }

    private void createQueryJsonFile() throws IOException {
        String queryJsonFileContent = buildQueryJson(WORKSPACE.getAbsolutePath(),
                TEST_CHROMOSOME_NAME, TEMP_CALLSET_JSON_FILE.getAbsolutePath(), TEMP_VID_JSON_FILE.getAbsolutePath());
        FileOutputStream outputStream = new FileOutputStream(TEMP_QUERY_JSON_FILE);
        byte[] strToBytes = queryJsonFileContent.getBytes();
        outputStream.write(strToBytes);
        outputStream.close();
    }
}
