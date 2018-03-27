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

package com.intel.genomicsdb.importer;

import com.intel.genomicsdb.GenomicsDBTestUtils;
import com.intel.genomicsdb.importer.extensions.CallSetMapExtensions;
import com.intel.genomicsdb.importer.model.ChromosomeInterval;
import com.intel.genomicsdb.model.*;
import com.intel.genomicsdb.reader.GenomicsDBFeatureReader;
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
import org.testng.annotations.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.Optional;


public final class GenomicsDBImporterSpec implements CallSetMapExtensions {
    private static final String TEST_CHROMOSOME_NAME = "1";
    private static final File WORKSPACE = new File("__workspace");
    private File tempVidJsonFile;
    private File tempCallsetJsonFile;
    private File tempVcfHeaderFile;

    @Test(testName = "genomicsdb importer with an interval and multiple GVCFs",
            dataProvider = "vcfFiles",
            dataProviderClass = GenomicsDBTestUtils.class)
    public void testMultiGVCFInputs(Map<String, FeatureReader<VariantContext>> sampleToReaderMap)
            throws IOException, InterruptedException {
        GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB =
                this.generateSortedCallSetMap(sampleToReaderMap, true, false);

        GenomicsDBImporter importer = getGenomicsDBImporterForMultipleImport();

        importer.executeImport();

        importer.writeCallsetMapJSONFile(tempCallsetJsonFile.getAbsolutePath(),
                callsetMappingPB);
        Assert.assertEquals(importer.isDone(), true);
        File callsetFile = new File(tempCallsetJsonFile.getAbsolutePath());
        Assert.assertTrue(callsetFile.exists());
    }

    @Test(testName = "genomicsdb importer outputs merged headers as a JSON file",
            dataProvider = "vcfFiles",
            dataProviderClass = GenomicsDBTestUtils.class)
    public void testVidMapJSONOutput(Map<String, FeatureReader<VariantContext>> sampleToReaderMap)
            throws IOException, InterruptedException {
        GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB_A =
                this.generateSortedCallSetMap(sampleToReaderMap, true, false);

        Set<VCFHeaderLine> mergedHeader = createMergedHeader(sampleToReaderMap);

        GenomicsDBImporter importer = getGenomicsDBImporterForMultipleImport();

        importer.writeVidMapJSONFile(tempVidJsonFile.getAbsolutePath(),
                importer.generateVidMapFromMergedHeader(mergedHeader));
        importer.writeVcfHeaderFile(tempVcfHeaderFile.getAbsolutePath(), mergedHeader);
        importer.writeCallsetMapJSONFile(tempCallsetJsonFile.getAbsolutePath(),
                callsetMappingPB_A);
        importer.executeImport();

        Assert.assertEquals(importer.isDone(), true);
        Assert.assertEquals(tempVidJsonFile.isFile(), true);
        Assert.assertEquals(tempCallsetJsonFile.isFile(), true);
        Assert.assertEquals(tempVcfHeaderFile.isFile(), true);

        JSONParser parser = new JSONParser();

        try {
            FileReader fileReader = new FileReader(tempCallsetJsonFile.getAbsolutePath());
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

        this.generateSortedCallSetMap(sampleToReaderMap, true, false);
    }

    @Test(testName = "should run parallel import with multiple chromosome intervals and one GVCF")
    public void testShouldRunParallelImportWithMultipleChromosomeIntervalsAndOneGvcf() throws IOException, InterruptedException {
        //Given
        GenomicsDBImporter importer = getGenomicsDBImporterForMultipleImport();

        //When
        importer.executeImport();

        //Then
        Assert.assertEquals(tempVidJsonFile.isFile(), true);
        Assert.assertEquals(tempCallsetJsonFile.isFile(), true);

        String referenceGenome = "tests/inputs/chr1_10MB.fasta.gz";

        Coordinates.ContigInterval interval0 = Coordinates.ContigInterval.newBuilder().setBegin(0).setEnd(50000).setContig("").build();
        Coordinates.GenomicsDBColumnInterval columnInterval = Coordinates.GenomicsDBColumnInterval.newBuilder()
                .setContigInterval(interval0).build();
        Coordinates.GenomicsDBColumnOrInterval columnRange = Coordinates.GenomicsDBColumnOrInterval.newBuilder()
                .setColumnInterval(columnInterval).build();
        List<Coordinates.GenomicsDBColumnOrInterval> columnRangeLists = new ArrayList<>(2);
        columnRangeLists.add(columnRange);
        GenomicsDBExportConfiguration.GenomicsDBColumnOrIntervalList columnRanges =
                GenomicsDBExportConfiguration.GenomicsDBColumnOrIntervalList.newBuilder()
                        .addAllColumnOrIntervalList(columnRangeLists).build();
        GenomicsDBExportConfiguration.RowRangeList.Builder rowRange = GenomicsDBExportConfiguration.RowRangeList.newBuilder()
                .addRangeList(GenomicsDBExportConfiguration.RowRange.newBuilder().setHigh(3).setLow(0));

        GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration = GenomicsDBExportConfiguration.ExportConfiguration.newBuilder()
                .setWorkspace(WORKSPACE.getAbsolutePath()).setReferenceGenome(referenceGenome)
                .setVidMappingFile(tempVidJsonFile.getAbsolutePath())
                .setCallsetMappingFile(tempCallsetJsonFile.getAbsolutePath()).setProduceGTField(true)
                .setScanFull(true)
                .addQueryColumnRanges(columnRanges)
                .addQueryRowRanges(rowRange)
                .addAttributes("GT")
                .build();

        GenomicsDBFeatureReader reader = new GenomicsDBFeatureReader<>(exportConfiguration, new BCF2Codec(), Optional.empty());

        CloseableTribbleIterator<VariantContext> gdbIterator = reader.iterator();
        List<VariantContext> varCtxList = new ArrayList<>();
        while (gdbIterator.hasNext()) varCtxList.add(gdbIterator.next());

        assert(varCtxList.size() == 2);
        assert(varCtxList.get(0).getStart() == 12141);
        assert(varCtxList.get(1).getStart() == 17385);
    }

    @Test(testName = "should be able to query using specific array name")
    public void testShouldBeAbleToQueryUsingSpecificArrayName() throws IOException, InterruptedException {
        //Given
        GenomicsDBImporter importer = getGenomicsDBImporterForMultipleImport();
        importer.executeImport();
        Assert.assertEquals(tempVidJsonFile.isFile(), true);
        Assert.assertEquals(tempCallsetJsonFile.isFile(), true);
        String referenceGenome = "tests/inputs/chr1_10MB.fasta.gz";

        //When
        GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration = GenomicsDBExportConfiguration.ExportConfiguration.newBuilder()
                .setWorkspace(WORKSPACE.getAbsolutePath())
                .setReferenceGenome(referenceGenome)
                .setVidMappingFile(tempVidJsonFile.getAbsolutePath())
                .setCallsetMappingFile(tempCallsetJsonFile.getAbsolutePath()).setProduceGTField(true)
                .setScanFull(true)
                .setArray("1#17000#18000")
                .build();

        GenomicsDBFeatureReader reader = new GenomicsDBFeatureReader<>(exportConfiguration, new BCF2Codec(), Optional.empty());

        CloseableTribbleIterator<VariantContext> gdbIterator = reader.iterator();
        List<VariantContext> varCtxList = new ArrayList<>();
        while (gdbIterator.hasNext()) varCtxList.add(gdbIterator.next());

        //Then
        assert(varCtxList.size() == 1);
        assert(varCtxList.get(0).getStart() == 17385);
    }

    @Test(testName = "should be able to query using specific chromosome interval")
    public void testShouldBeAbleToQueryUsingSpecificChromosomeInterval() throws IOException, InterruptedException {
        //Given
        GenomicsDBImporter importer = getGenomicsDBImporterForMultipleImport();
        importer.executeImport();
        Assert.assertEquals(tempVidJsonFile.isFile(), true);
        Assert.assertEquals(tempCallsetJsonFile.isFile(), true);

        //When
        String referenceGenome = "tests/inputs/chr1_10MB.fasta.gz";

        GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration = GenomicsDBExportConfiguration.ExportConfiguration.newBuilder()
                .setWorkspace(WORKSPACE.getAbsolutePath())
                .setReferenceGenome(referenceGenome)
                .setVidMappingFile(tempVidJsonFile.getAbsolutePath())
                .setCallsetMappingFile(tempCallsetJsonFile.getAbsolutePath())
                .setProduceGTField(true)
                .setScanFull(true)
                .build();

        GenomicsDBFeatureReader reader = new GenomicsDBFeatureReader<>(exportConfiguration, new BCF2Codec(), Optional.empty());

        CloseableTribbleIterator<VariantContext> gdbIterator = reader.query("1", 17000, 18000);
        List<VariantContext> varCtxList = new ArrayList<>();
        while (gdbIterator.hasNext()) varCtxList.add(gdbIterator.next());

        //Then
        assert(varCtxList.size() == 1);
        assert(varCtxList.get(0).getStart() == 17385);
    }

    @BeforeMethod
    public void cleanUpBefore() throws IOException {
        tempVidJsonFile = File.createTempFile("generated_vidmap", ".json"); //new File("generated_vidmap.json");
        tempCallsetJsonFile = File.createTempFile("generated_callsetmap", ".json"); //new File("generated_callsetmap.json");
        tempVcfHeaderFile = File.createTempFile("vcfheader", ".vcf"); //new File("vcfheader.vcf");
    }

    @AfterMethod
    public void cleanUpAfter() throws IOException {
        if (tempVidJsonFile.exists()) FileUtils.deleteQuietly(tempVidJsonFile);
        if (tempCallsetJsonFile.exists()) FileUtils.deleteQuietly(tempCallsetJsonFile);
        if (WORKSPACE.exists()) FileUtils.deleteDirectory(WORKSPACE);
    }

    private Set<VCFHeaderLine> createMergedHeader(Map<String, FeatureReader<VariantContext>> variantReaders) {
        List<VCFHeader> headers = new ArrayList<>();
        for (Map.Entry<String, FeatureReader<VariantContext>> variant : variantReaders.entrySet()) {
            headers.add((VCFHeader) variant.getValue().getHeader());
        }
        return VCFUtils.smartMergeHeaders(headers, true);
    }

    private GenomicsDBImporter getGenomicsDBImporterForMultipleImport() throws FileNotFoundException {
        String[] args = ("-L 1:12000-13000 -L 1:17000-18000 " +
                "-w " + WORKSPACE.getAbsolutePath() + " " +
                "--size_per_column_partition 16384 " +
                "--segment_size 10485760 " +
                "--vidmap-output " + tempVidJsonFile.getAbsolutePath() + " " +
                "--callset-output " + tempCallsetJsonFile.getAbsolutePath() + " " +
                "tests/inputs/vcfs/t0.vcf.gz").split(" ");
        CommandLineImportConfig config = new CommandLineImportConfig("TestGenomicsDBImporterWithMergedVCFHeader", args);
        return new GenomicsDBImporter(config);
    }
}
