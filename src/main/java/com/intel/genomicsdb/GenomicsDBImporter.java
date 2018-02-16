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

import com.googlecode.protobuf.format.JsonFormat;
import com.intel.genomicsdb.model.BaseImportConfig;
import com.intel.genomicsdb.model.VCFHeaderToFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

import static com.googlecode.protobuf.format.JsonFormat.printToString;

/**
 * Java wrapper for vcf2tiledb - imports VCFs into TileDB/GenomicsDB.
 * All vid information is assumed to be set correctly by the user (JSON files)
 */
public class GenomicsDBImporter {
    private static final String DEFAULT_ARRAY_NAME = "genomicsdb_array";
    private static final int DEFAULT_TILEDB_CELLS_PER_TILE = 1000;
    //Allele specific annotation fields
    private static final HashSet<String> mRLengthHistogramFieldsWithFloatBins = new HashSet<>(Arrays.asList(
            "AS_RAW_BaseQRankSum",
            "AS_RAW_MQRankSum",
            "AS_RAW_ReadPosRankSum"
    ));
    private static final HashSet<String> mRLengthTwoDFloatVectorFields = new HashSet<>(Collections.singletonList(
            "AS_RAW_MQ"
    ));
    private static final HashSet<String> mRLengthTwoDIntVectorFields = new HashSet<>(Collections.singletonList(
            "AS_SB_TABLE"
    ));
    static long mDefaultBufferCapacity = 20480; //20KB

    static {
        try {
            boolean loaded = GenomicsDBUtils.loadLibrary();
            if (!loaded) throw new GenomicsDBException("Could not load genomicsdb native library");
        } catch (UnsatisfiedLinkError ule) {
            throw new GenomicsDBException("Could not load genomicsdb native library", ule);
        }
    }

  /*
   * JNI functions
   */

    private String mLoaderJSONFile = null;
    private int mRank = 0;
    private long mLbRowIdx = 0;
    private long mUbRowIdx = Long.MAX_VALUE - 1;
    //For buffered streams
    private boolean mContainsBufferStreams = false;
    private long mGenomicsDBImporterObjectHandle = 0;
    private ArrayList<GenomicsDBImporterStreamWrapper> mBufferStreamWrapperVector = null;
    private boolean mIsLoaderSetupDone = false;
    private long[] mExhaustedBufferStreamIdentifiers = null;
    private long mNumExhaustedBufferStreams = 0;
    //Done flag - useful only for buffered streams
    private boolean mDone = false;
    //JSON object that specifies callset/sample name to row_idx mapping in the buffer
    private JSONObject mCallsetMappingJSON = null;
    private boolean mUsingVidMappingProtoBuf = false;

    /**
     * Default Constructor
     */
    public GenomicsDBImporter() {
    }

    /**
     * Constructor
     *
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     */
    public GenomicsDBImporter(String loaderJSONFile) {
        initialize(loaderJSONFile, 0, 0,
                Long.MAX_VALUE - 1);
    }

    /**
     * Constructor
     *
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank           Rank of this process (TileDB/GenomicsDB partition idx)
     */
    public GenomicsDBImporter(String loaderJSONFile, int rank) {
        initialize(loaderJSONFile, rank, 0, Long.MAX_VALUE - 1);
    }
    /**
     * Constructor
     *
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank           Rank of this process (TileDB/GenomicsDB partition idx)
     * @param lbRowIdx       Smallest row idx which should be imported by this object
     * @param ubRowIdx       Largest row idx which should be imported by this object
     */
    public GenomicsDBImporter(String loaderJSONFile, int rank, long lbRowIdx, long ubRowIdx) {
        initialize(loaderJSONFile, rank, lbRowIdx, ubRowIdx);
    }
    /**
     * Constructor to create required data structures from a list
     * of GVCF files and a chromosome interval. This constructor
     * is developed specifically for GATK4 GenomicsDBImport tool.
     *
     * @param sampleToVCMap          Variant Readers objects of the input GVCF files
     * @param mergedHeader           Headers from all input GVCF files merged into one
     * @param chromosomeInterval     Chromosome interval to traverse input VCFs
     * @param workspace              TileDB workspace
     * @param arrayname              TileDB array name
     * @param sizePerColumnPartition sizePerColumnPartition in bytes
     * @param segmentSize            segmentSize in bytes
     * @throws IOException FeatureReader.query can throw an IOException if invalid file is used
     */
    public GenomicsDBImporter(Map<String, FeatureReader<VariantContext>> sampleToVCMap,
                              Set<VCFHeaderLine> mergedHeader,
                              ChromosomeInterval chromosomeInterval,
                              String workspace,
                              String arrayname,
                              Long sizePerColumnPartition,
                              Long segmentSize) throws IOException {
        this(sampleToVCMap, mergedHeader, chromosomeInterval,
                workspace, arrayname, sizePerColumnPartition, segmentSize,
                (long) 0, Long.MAX_VALUE - 1, true);
    }
    /**
     * Constructor with an GenomicsDB import configuration protocol buffer
     * structure. Avoids passing long list of parameters. This constructor
     * is developed specifically for GATK4 GenomicsDBImport tool.
     *
     * @param sampleToVCMap             Variant Readers objects of the input GVCF files
     * @param mergedHeader              Set of VCFHeaderLine from the merged header across all input files
     * @param chromosomeInterval        Chromosome interval to traverse input VCFs
     * @param importConfiguration       Protobuf configuration object containing related input
     *                                  parameters, filenames, etc.
     * @param validateSampleToReaderMap Check validity of sampleToreaderMap entries
     * @throws IOException Throws file IO exception.
     */
    public GenomicsDBImporter(Map<String, FeatureReader<VariantContext>> sampleToVCMap,
                              Set<VCFHeaderLine> mergedHeader,
                              ChromosomeInterval chromosomeInterval,
                              boolean validateSampleToReaderMap,
                              GenomicsDBImportConfiguration.ImportConfiguration importConfiguration)
            throws IOException {
        this(sampleToVCMap,
                mergedHeader,
                chromosomeInterval,
                importConfiguration.getColumnPartitions(0).getWorkspace(),
                importConfiguration.getColumnPartitions(0).getArray(),
                importConfiguration.getSizePerColumnPartition(),
                importConfiguration.getSegmentSize(),
                importConfiguration.getGatk4IntegrationParameters().getLowerSampleIndex(),
                importConfiguration.getGatk4IntegrationParameters().getUpperSampleIndex(),
                false,  //useSamplesInOrderProvided
                importConfiguration.getFailIfUpdating(),
                validateSampleToReaderMap
        );
    }

    /**
     * Constructor to create required data structures from a list
     * of GVCF files and a chromosome interval. This constructor
     * is developed specifically for GATK4 GenomicsDBImport tool.
     *
     * @param sampleToReaderMap         Variant Readers objects of the input GVCF files
     * @param mergedHeader              Headers from all input GVCF files merged into one
     * @param chromosomeInterval        Chromosome interval to traverse input VCFs
     * @param workspace                 TileDB workspace
     * @param arrayname                 TileDB array name
     * @param sizePerColumnPartition    sizePerColumnPartition in bytes
     * @param segmentSize               segmentSize in bytes
     * @param lbRowIdx                  Smallest row idx which should be imported by this object
     * @param ubRowIdx                  Largest row idx which should be imported by this object
     * @param validateSampleToReaderMap Check validity of sampleToreaderMap entries
     * @throws IOException FeatureReader.query can throw an IOException if invalid file is used
     */
    public GenomicsDBImporter(Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
                              Set<VCFHeaderLine> mergedHeader,
                              ChromosomeInterval chromosomeInterval,
                              String workspace,
                              String arrayname,
                              Long sizePerColumnPartition,
                              Long segmentSize,
                              Long lbRowIdx,
                              Long ubRowIdx,
                              boolean validateSampleToReaderMap) throws IOException {
        this(sampleToReaderMap, mergedHeader, chromosomeInterval,
                workspace, arrayname, sizePerColumnPartition, segmentSize, lbRowIdx, ubRowIdx,
                false, validateSampleToReaderMap);
    }
    /**
     * Constructor to create required data structures from a list
     * of GVCF files and a chromosome interval. This constructor
     * is developed specifically for GATK4 GenomicsDBImport tool.
     *
     * @param sampleToReaderMap         Variant Readers objects of the input GVCF files
     * @param mergedHeader              Headers from all input GVCF files merged into one
     * @param chromosomeInterval        Chromosome interval to traverse input VCFs
     * @param workspace                 TileDB workspace
     * @param arrayname                 TileDB array name
     * @param sizePerColumnPartition    sizePerColumnPartition in bytes
     * @param segmentSize               segmentSize in bytes
     * @param lbRowIdx                  Smallest row idx which should be imported by this object
     * @param ubRowIdx                  Largest row idx which should be imported by this object
     * @param useSamplesInOrderProvided if true, don't sort samples, instead
     *                                  use in the the order provided
     * @param validateSampleToReaderMap Check validity of sampleToreaderMap entries
     * @throws IOException FeatureReader.query can throw an IOException if invalid file is used
     */
    public GenomicsDBImporter(Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
                              Set<VCFHeaderLine> mergedHeader,
                              ChromosomeInterval chromosomeInterval,
                              String workspace,
                              String arrayname,
                              Long sizePerColumnPartition,
                              Long segmentSize,
                              Long lbRowIdx,
                              Long ubRowIdx,
                              boolean useSamplesInOrderProvided,
                              boolean validateSampleToReaderMap) throws IOException {
        this(sampleToReaderMap, mergedHeader, chromosomeInterval,
                workspace, arrayname, sizePerColumnPartition, segmentSize,
                lbRowIdx, ubRowIdx, useSamplesInOrderProvided, false, validateSampleToReaderMap);
    }
    /**
     * Constructor to create required data structures from a list
     * of GVCF files and a chromosome interval. This constructor
     * is developed specifically for GATK4 GenomicsDBImport tool.
     *
     * @param sampleToReaderMap         Variant Readers objects of the input GVCF files
     * @param mergedHeader              Set of VCFHeaderLine from the merged header across all input files
     * @param chromosomeInterval        Chromosome interval to traverse input VCFs
     * @param workspace                 TileDB workspace
     * @param arrayname                 TileDB array name
     * @param sizePerColumnPartition    sizePerColumnPartition in bytes
     * @param segmentSize               segmentSize in bytes
     * @param lbRowIdx                  Smallest row idx which should be imported by this object
     * @param ubRowIdx                  Largest row idx which should be imported by this object
     * @param useSamplesInOrderProvided if true, don't sort samples, instead use in the the order
     *                                  provided
     * @param failIfUpdating            if true, fail if updating an existing array
     * @param validateSampleToReaderMap Check validity of sampleToreaderMap entries
     * @throws IOException when load into TileDB array fails
     */
    public GenomicsDBImporter(Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
                              Set<VCFHeaderLine> mergedHeader,
                              ChromosomeInterval chromosomeInterval,
                              String workspace,
                              String arrayname,
                              Long sizePerColumnPartition,
                              Long segmentSize,
                              Long lbRowIdx,
                              Long ubRowIdx,
                              boolean useSamplesInOrderProvided,
                              boolean failIfUpdating,
                              boolean validateSampleToReaderMap) throws IOException {
        this(sampleToReaderMap, mergedHeader, chromosomeInterval,
                workspace, arrayname, sizePerColumnPartition, segmentSize,
                lbRowIdx, ubRowIdx, useSamplesInOrderProvided, failIfUpdating, 0, validateSampleToReaderMap);
    }
    /**
     * Constructor to create required data structures from a list
     * of GVCF files and a chromosome interval. This constructor
     * is developed specifically for GATK4 GenomicsDBImport tool.
     *
     * @param sampleToReaderMap               Feature Readers objects corresponding to input GVCF files
     * @param mergedHeader                    Set of VCFHeaderLine from the merged header across all input files
     * @param chromosomeInterval              Chromosome interval to traverse input VCFs
     * @param workspace                       TileDB workspace
     * @param arrayname                       TileDB array name
     * @param vcfBufferSizePerColumnPartition vcfBufferSizePerColumnPartition in bytes
     * @param segmentSize                     segmentSize in bytes
     * @param lbRowIdx                        Smallest row idx which should be imported by this object
     * @param ubRowIdx                        Largest row idx which should be imported by this object
     * @param useSamplesInOrderProvided       if true, don't sort samples, instead use in the the order
     *                                        provided
     * @param failIfUpdating                  if true, fail if updating an existing array
     * @param rank                            Rank of object - corresponds to the partition index in the loader
     * @param validateSampleToReaderMap       Check validity of sampleToreaderMap entries
     * @throws IOException when load into TileDB array fails
     */
    public GenomicsDBImporter(Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
                              Set<VCFHeaderLine> mergedHeader,
                              ChromosomeInterval chromosomeInterval,
                              String workspace,
                              String arrayname,
                              Long vcfBufferSizePerColumnPartition,
                              Long segmentSize,
                              Long lbRowIdx,
                              Long ubRowIdx,
                              boolean useSamplesInOrderProvided,
                              boolean failIfUpdating,
                              int rank,
                              boolean validateSampleToReaderMap) throws IOException, IllegalArgumentException {
        this(sampleToReaderMap,
                mergedHeader,
                chromosomeInterval,
                workspace,
                arrayname,
                vcfBufferSizePerColumnPartition,
                segmentSize,
                lbRowIdx,
                ubRowIdx,
                useSamplesInOrderProvided,
                failIfUpdating,
                rank,
                validateSampleToReaderMap,
                true);
    }
    /**
     * Constructor to create required data structures from a list
     * of GVCF files and a chromosome interval. This constructor
     * is developed specifically for GATK4 GenomicsDBImport tool.
     *
     * @param sampleToReaderMap               Feature Readers objects corresponding to input GVCF files
     * @param mergedHeader                    Set of VCFHeaderLine from the merged header across all input files
     * @param chromosomeInterval              Chromosome interval to traverse input VCFs
     * @param workspace                       TileDB workspace
     * @param arrayname                       TileDB array name
     * @param vcfBufferSizePerColumnPartition vcfBufferSizePerColumnPartition in bytes
     * @param segmentSize                     segmentSize in bytes
     * @param lbRowIdx                        Smallest row idx which should be imported by this object
     * @param ubRowIdx                        Largest row idx which should be imported by this object
     * @param useSamplesInOrderProvided       if true, don't sort samples, instead use in the the order
     *                                        provided
     * @param failIfUpdating                  if true, fail if updating an existing array
     * @param rank                            Rank of object - corresponds to the partition index in the loader
     * @param validateSampleToReaderMap       Check validity of sampleToreaderMap entries
     * @param passAsVcf                       Use the VCF format to pass data from Java to C++
     * @throws IOException when load into TileDB array fails
     */
    public GenomicsDBImporter(Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
                              Set<VCFHeaderLine> mergedHeader,
                              ChromosomeInterval chromosomeInterval,
                              String workspace,
                              String arrayname,
                              Long vcfBufferSizePerColumnPartition,
                              Long segmentSize,
                              Long lbRowIdx,
                              Long ubRowIdx,
                              boolean useSamplesInOrderProvided,
                              boolean failIfUpdating,
                              int rank,
                              boolean validateSampleToReaderMap,
                              boolean passAsVcf) throws IOException, IllegalArgumentException {
        // Mark this flag so that protocol buffer based vid
        // and callset map are propagated to C++ GenomicsDBImporter
        mUsingVidMappingProtoBuf = true;

        GenomicsDBImportConfiguration.ImportConfiguration importConfiguration =
                createImportConfiguration(workspace, arrayname, vcfBufferSizePerColumnPartition,
                        segmentSize, failIfUpdating);


        GenomicsDBVidMapProto.VidMappingPB vidMapPB = generateVidMapFromMergedHeader(mergedHeader);

        GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMapPB =
                generateSortedCallSetMap(sampleToReaderMap, useSamplesInOrderProvided,
                        validateSampleToReaderMap, lbRowIdx);

        File importJSONFile = dumpTemporaryLoaderJSONFile(importConfiguration, "");

        initialize(importJSONFile.getAbsolutePath(), rank, lbRowIdx, ubRowIdx);

        mGenomicsDBImporterObjectHandle =
                jniInitializeGenomicsDBImporterObject(mLoaderJSONFile, mRank, lbRowIdx, ubRowIdx);

        jniCopyVidMap(mGenomicsDBImporterObjectHandle, vidMapPB.toByteArray());
        jniCopyCallsetMap(mGenomicsDBImporterObjectHandle, callsetMapPB.toByteArray());

        for (GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap sampleToIDMap :
                callsetMapPB.getCallsetsList()) {

            String sampleName = sampleToIDMap.getSampleName();

            FeatureReader<VariantContext> featureReader = sampleToReaderMap.get(sampleName);

            CloseableIterator<VariantContext> iterator = featureReader.query(chromosomeInterval.getContig(),
                    chromosomeInterval.getStart(), chromosomeInterval.getEnd());

            addSortedVariantContextIterator(
                    sampleToIDMap.getStreamName(),
                    (VCFHeader) featureReader.getHeader(),
                    iterator,
                    importConfiguration.getSizePerColumnPartition(),
                    passAsVcf ? VariantContextWriterBuilder.OutputType.VCF_STREAM
                            : VariantContextWriterBuilder.OutputType.BCF_STREAM,
                    null);
        }
    }

    /**
     * Obtain the chromosome intervals for the column partition specified in the loader JSON file
     * identified by the rank. The information is returned as a string in JSON format
     * {
     * "contigs": [
     * { "chr1": [ 100, 200] },
     * { "chr2": [ 500, 600] }
     * ]
     * }
     *
     * @param loaderJSONFile path to loader JSON file
     * @param rank           rank/partition index
     * @return chromosome intervals for the queried column partition in JSON format
     */
    private static native String jniGetChromosomeIntervalsForColumnPartition(
            final String loaderJSONFile,
            final int rank);

    /**
     * Create TileDB workspace
     *
     * @param workspace path to workspace directory
     * @return status 0 = workspace created,
     * -1 = path was not a directory,
     * -2 = failed to create workspace,
     * 1 = existing directory, nothing changed
     */
    private static native int jniCreateTileDBWorkspace(final String workspace);

    /**
     * Consolidate TileDB array
     *
     * @param workspace path to workspace directory
     * @param arrayName array name
     */
    private static native void jniConsolidateTileDBArray(final String workspace, final String arrayName);

    private static VCFHeaderToFile resolveVCFHeaderToFile(final BaseImportConfig baseImportConfig) {
        List<VCFHeader> headers = new ArrayList<>();
        ArrayList<String> sampleNames = new ArrayList<>();
        Map<String, String> sampleNameToFileName = new LinkedHashMap<>();

        //Get merged header first
        for (String file : baseImportConfig.getFiles()) {
            AbstractFeatureReader<VariantContext, LineIterator> reader =
                    AbstractFeatureReader.getFeatureReader(file, new VCFCodec(), false);
            headers.add((VCFHeader) reader.getHeader());
            final String sampleName = ((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0);
            sampleNames.add(sampleName);
            sampleNameToFileName.put(sampleName, file);
            //Hopefully, GC kicks in and frees resources assigned to reader
        }

        return new VCFHeaderToFile(headers, sampleNames, sampleNameToFileName);
    }

    /**
     * Import multiple chromosome interval
     * @param baseImportConfig
     */
    public static void parallelImport(final BaseImportConfig baseImportConfig) throws IOException {
        VCFHeaderToFile vcfHeaderToFile = resolveVCFHeaderToFile(baseImportConfig);

        //Merge headers
        Set<VCFHeaderLine> mergedHeader = VCFUtils.smartMergeHeaders(vcfHeaderToFile.getHeaders(), true);

        //This sorts the list sampleNames if !useSamplesInOrder
        //Why you should use this? If you are writing multiple partitions in different machines,
        //you must have consistent ordering of samples across partitions. If file order is different
        //in different processes, then set useSamplesInOrder to false and let the sort in
        //generateSortedCallSetMap ensure consistent ordering across samples
        GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB =
                GenomicsDBImporter.generateSortedCallSetMap(vcfHeaderToFile.getSampleNames(), baseImportConfig.isUseSamplesInOrderProvided());

        //Write out callset map if needed
        if (baseImportConfig.getCallsetOutputFilepath() != null)
            GenomicsDBImporter.writeCallsetMapJSONFile(baseImportConfig.getCallsetOutputFilepath(), callsetMappingPB);

        //Write out vidmap if needed
        if (baseImportConfig.getVidmapOutputFilepath() != null)
            GenomicsDBImporter.writeVidMapJSONFile(baseImportConfig.getVidmapOutputFilepath(), mergedHeader);

        //write out merged header if needed
        if (baseImportConfig.getVcfHeaderOutputFilepath() != null)
            GenomicsDBImporter.writeVcfHeaderFile(baseImportConfig.getVcfHeaderOutputFilepath(), mergedHeader);

        //Iterate over sorted sample list in batches
        //TODO: make sample list iteration in parallel
        for (int i = 0; i < vcfHeaderToFile.getSampleNames().size(); i += baseImportConfig.getBatchSize()) {
            final int index = i;
            final Map<String, FeatureReader<VariantContext>> sampleToReaderMap =
                    createSampleToReaderMap(baseImportConfig, vcfHeaderToFile, index);

            List<Boolean> result = baseImportConfig.getChromosomeIntervalList().parallelStream().map(chromosomeInterval -> {
                try {
                    GenomicsDBImporter importer = createImporter(
                            baseImportConfig, vcfHeaderToFile, mergedHeader, index, sampleToReaderMap, chromosomeInterval);
                    return importer.importBatch();
                } catch (Exception ex) {
                    throw new IllegalStateException("There was an unhandled exception during chromosome interval import.");
                }
            }).collect(Collectors.toList());

            if(result.contains(false))
               throw new IllegalStateException("There was an unhandled exception during chromosome interval import.");
        }
    }

    private static Map<String, FeatureReader<VariantContext>> createSampleToReaderMap(BaseImportConfig baseImportConfig, VCFHeaderToFile vcfHeaderToFile, int index) {
        final Map<String, FeatureReader<VariantContext>> sampleToReaderMap = new LinkedHashMap<>();
        for (int j = index; j < vcfHeaderToFile.getSampleNames().size() && j < index + baseImportConfig.getBatchSize(); ++j) {
            final String sampleName = vcfHeaderToFile.getSampleNames().get(j);
            assert vcfHeaderToFile.getSampleNameToFileName().containsKey(sampleName);
            AbstractFeatureReader<VariantContext, LineIterator> reader =
                    AbstractFeatureReader.getFeatureReader(vcfHeaderToFile.getSampleNameToFileName().get(sampleName), new VCFCodec(), false);
            assert sampleName.equals(((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0));
            sampleToReaderMap.put(sampleName, reader);
        }
        return sampleToReaderMap;
    }

    private static synchronized GenomicsDBImporter createImporter(final BaseImportConfig baseImportConfig, final VCFHeaderToFile vcfHeaderToFile,
                                                                  final Set<VCFHeaderLine> mergedHeader, final int index,
                                                                  final Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
                                                                  final ChromosomeInterval chromInterval) throws IOException {
        return new GenomicsDBImporter(
                sampleToReaderMap, mergedHeader,
                new ChromosomeInterval(chromInterval.getContig(), chromInterval.getStart(), chromInterval.getEnd()),
                baseImportConfig.getWorkspace(), chromInterval.getContig(),
                baseImportConfig.getVcfBufferSizePerColumnPartition() * vcfHeaderToFile.getSampleNames().size(),
                baseImportConfig.getSegmentSize(), (long) index, (long) (index + baseImportConfig.getBatchSize() - 1),
                baseImportConfig.isUseSamplesInOrderProvided(), baseImportConfig.isFailIfUpdating(),
                baseImportConfig.getRank(), baseImportConfig.isValidateSampleToReaderMap(), baseImportConfig.isPassAsVcf());
    }

    /**
     * Create TileDB workspace
     *
     * @param workspace path to workspace directory
     * @return status 0 = workspace created,
     * -1 = path was not a directory,
     * -2 = failed to create workspace,
     * 1 = existing directory, nothing changed
     */
    public static int createTileDBWorkspace(final String workspace) {
        return jniCreateTileDBWorkspace(workspace);
    }

    /**
     * Consolidate TileDB array
     *
     * @param workspace path to workspace directory
     * @param arrayName array name
     */
    public static void consolidateTileDBArray(final String workspace, final String arrayName) {
        jniConsolidateTileDBArray(workspace, arrayName);
    }

    /**
     * Writes a JSON file from a set of VCF header lines. This
     * method is not called implicitly by GenomicsDB constructor.
     * It needs to be called explicitly if the user wants these objects
     * to be written. Called explicitly from GATK-4 GenomicsDBImport tool
     *
     * @param outputVidMapJSONFilePath Full path of file to be written
     * @param headerLines              Set of header lines
     * @throws FileNotFoundException PrintWriter throws this exception when file not found
     */
    static void writeVidMapJSONFile(String outputVidMapJSONFilePath, Set<VCFHeaderLine> headerLines) throws FileNotFoundException {
        GenomicsDBVidMapProto.VidMappingPB vidMappingPB = generateVidMapFromMergedHeader(headerLines);
        String vidMapJSONString = printToString(vidMappingPB);
        File vidMapJSONFile = new File(outputVidMapJSONFilePath);

        PrintWriter out = new PrintWriter(vidMapJSONFile);
        out.println(vidMapJSONString);
        out.close();
    }

    /**
     * Writes a VCF Header file from a set of VCF header lines. This
     * method is not called implicitly by GenomicsDB constructor.
     * It needs to be called explicitly if the user wants these objects
     * to be written. Called explicitly from GATK-4 GenomicsDBImport tool
     *
     * @param outputVcfHeaderFilePath Full path of file to be written
     * @param headerLines             Set of header lines
     * @throws FileNotFoundException PrintWriter throws this exception when file not found
     */
    static void writeVcfHeaderFile(String outputVcfHeaderFilePath, Set<VCFHeaderLine> headerLines) throws FileNotFoundException {
        File vcfHeaderFile = new File(outputVcfHeaderFilePath);
        final VCFHeader vcfHeader = new VCFHeader(headerLines);
        VariantContextWriter vcfWriter = new VariantContextWriterBuilder()
                .clearOptions()
                .setOutputFile(vcfHeaderFile)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .build();
        vcfWriter.writeHeader(vcfHeader);
        vcfWriter.close();
    }

    /**
     * Writes a JSON file from a vidmap protobuf object. This
     * method is not called implicitly by GenomicsDB constructor.
     * It needs to be called explicitly if the user wants these objects
     * to be written. Called explicitly from GATK-4 GenomicsDBImport tool
     *
     * @param outputCallsetMapJSONFilePath Full path of file to be written
     * @param callsetMappingPB             Protobuf callset map object
     * @throws FileNotFoundException PrintWriter throws this exception when file not found
     */
    static void writeCallsetMapJSONFile(String outputCallsetMapJSONFilePath,
                                               GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB) throws FileNotFoundException {
        String callsetMapJSONString = printToString(callsetMappingPB);
        File callsetMapJSONFile = new File(outputCallsetMapJSONFilePath);

        PrintWriter out = new PrintWriter(callsetMapJSONFile);
        out.println(callsetMapJSONString);
        out.close();
    }

    /**
     * Create a JSON file from the import configuration
     *
     * @param importConfiguration The configuration object
     * @param filename            File to dump the loader JSON to
     * @return New file (with the specified name) written to local storage
     * @throws FileNotFoundException PrintWriter throws this exception when file not found
     */
    private static File dumpTemporaryLoaderJSONFile(
            GenomicsDBImportConfiguration.ImportConfiguration importConfiguration,
            String filename) throws IOException {
        String loaderJSONString = JsonFormat.printToString(importConfiguration);

        File tempLoaderJSONFile = (filename.isEmpty()) ?
                File.createTempFile("loader_", ".json") :
                new File(filename);

        if (filename.isEmpty())
            tempLoaderJSONFile.deleteOnExit();

        PrintWriter out = new PrintWriter(tempLoaderJSONFile);
        out.println(loaderJSONString);
        out.close();
        return tempLoaderJSONFile;
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes.
     * <p>
     * Assume one sample per input GVCF file
     *
     * @param sampleToReaderMap         Variant Readers objects of the input GVCF files
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @param validateSampleToReaderMap
     * @return Mappings of callset (sample) names to TileDB rows
     */
    static GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(
            final Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
            boolean validateSampleToReaderMap,
            boolean useSamplesInOrderProvided) {
        return GenomicsDBImporter.generateSortedCallSetMap(sampleToReaderMap,
                useSamplesInOrderProvided, validateSampleToReaderMap, 0L);
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes.
     * <p>
     * Assume one sample per input GVCF file
     *
     * @param sampleToReaderMap         Variant Readers objects of the input GVCF files
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @param validateSampleMap         Check i) whether sample names are consistent
     *                                  with headers and ii) feature readers are valid
     *                                  in sampleToReaderMap
     * @param lbRowIdx                  Smallest row idx which should be imported by this object
     * @return Mappings of callset (sample) names to TileDB rows
     */
    private static GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(
            final Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
            boolean useSamplesInOrderProvided,
            boolean validateSampleMap,
            final long lbRowIdx) {
        if (!validateSampleMap) {
            return GenomicsDBImporter.generateSortedCallSetMap(
                    new ArrayList<>(sampleToReaderMap.keySet()), useSamplesInOrderProvided, lbRowIdx);
        }

        List<String> listOfSampleNames = new ArrayList<>(sampleToReaderMap.size());
        for (Map.Entry<String, FeatureReader<VariantContext>> mapObject :
                sampleToReaderMap.entrySet()) {

            String sampleName = mapObject.getKey();
            FeatureReader<VariantContext> featureReader = mapObject.getValue();

            if (featureReader == null) {
                throw new IllegalArgumentException("Null FeatureReader found for sample: " + sampleName);
            }

            VCFHeader header = (VCFHeader) featureReader.getHeader();
            List<String> sampleNamesInHeader = header.getSampleNamesInOrder();
            if (sampleNamesInHeader.size() > 1) {
                StringBuilder messageBuilder = new StringBuilder("Multiple samples ");
                for (String name : sampleNamesInHeader) {
                    messageBuilder.append(name).append(" ");
                }
                messageBuilder.append(" appear in header for ").append(sampleName);
                throw new IllegalArgumentException(messageBuilder.toString());
            } else {
                if (!sampleName.equals(sampleNamesInHeader.get(0))) {
                    System.err.println("Sample " + header + " does not match with " +
                            sampleNamesInHeader.get(0) + " in header");
                }
            }
            listOfSampleNames.add(sampleName);
        }
        return GenomicsDBImporter.generateSortedCallSetMap(listOfSampleNames, useSamplesInOrderProvided,
                lbRowIdx);
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes.
     *
     * @param sampleNames list of sample names
     * @return Mappings of callset (sample) names to TileDB rows
     */
    public static GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(
            List<String> sampleNames) {
        return GenomicsDBImporter.generateSortedCallSetMap(
                sampleNames,
                false, 0L);
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes.
     *
     * @param sampleNames               list of sample names
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @return Mappings of callset (sample) names to TileDB rows
     */
    static GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(
            List<String> sampleNames,
            boolean useSamplesInOrderProvided) {
        return GenomicsDBImporter.generateSortedCallSetMap(
                sampleNames,
                useSamplesInOrderProvided, 0L);
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes. This method is synchronized
     * to block multiple invocations (if by any chance) disturb the order in
     * which TileDB row indexes are generated
     *
     * @param sampleNames               list of sample names
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @param lbRowIdx                  Smallest row idx which should be imported by this object
     * @return Mappings of callset (sample) names to TileDB rows
     */
    private static GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(
            List<String> sampleNames,
            boolean useSamplesInOrderProvided,
            final long lbRowIdx) {
        if (!useSamplesInOrderProvided)
            Collections.sort(sampleNames);

        GenomicsDBCallsetsMapProto.CallsetMappingPB.Builder callsetMapBuilder =
                GenomicsDBCallsetsMapProto.CallsetMappingPB.newBuilder();

        long tileDBRowIndex = lbRowIdx;

        for (String sampleName : sampleNames) {
            GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap.Builder idMapBuilder =
                    GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap.newBuilder();

            idMapBuilder
                    .setSampleName(sampleName)
                    .setRowIdx(tileDBRowIndex++)
                    .setIdxInFile(0)
                    .setStreamName(sampleName + "_stream");

            GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap sampleIDToTileDBIDMap =
                    idMapBuilder.build();

            callsetMapBuilder.addCallsets(sampleIDToTileDBIDMap);
        }
        return callsetMapBuilder.build();
    }

    /**
     * Generate the ProtoBuf data structure for vid mapping
     * Also remember, contigs are 1-based which means
     * TileDB column offsets should start from 1
     *
     * @param mergedHeader Header from all input GVCFs are merged to
     *                     create one vid map for all
     * @return a vid map containing all field names, lengths and types
     * from the merged GVCF header
     */
    private static GenomicsDBVidMapProto.VidMappingPB generateVidMapFromMergedHeader(
            Set<VCFHeaderLine> mergedHeader) {

        List<GenomicsDBVidMapProto.InfoField> infoFields = new ArrayList<>();
        List<GenomicsDBVidMapProto.Chromosome> contigs = new ArrayList<>();

        //ID field
        GenomicsDBVidMapProto.InfoField.Builder IDFieldBuilder =
                GenomicsDBVidMapProto.InfoField.newBuilder();
        GenomicsDBVidMapProto.FieldLengthDescriptorComponentPB.Builder lengthDescriptorComponentBuilder =
                GenomicsDBVidMapProto.FieldLengthDescriptorComponentPB.newBuilder();
        lengthDescriptorComponentBuilder.setVariableLengthDescriptor("var");
        IDFieldBuilder.setName("ID").addType("char").addLength(lengthDescriptorComponentBuilder.build());
        infoFields.add(IDFieldBuilder.build());

        int dpIndex = -1;
        long columnOffset = 0L;

        for (VCFHeaderLine headerLine : mergedHeader) {
            GenomicsDBVidMapProto.InfoField.Builder infoBuilder =
                    GenomicsDBVidMapProto.InfoField.newBuilder();

            if (headerLine instanceof VCFFormatHeaderLine) {
                VCFFormatHeaderLine formatHeaderLine = (VCFFormatHeaderLine) headerLine;
                boolean isGT = formatHeaderLine.getID().equals(VCFConstants.GENOTYPE_KEY);
                String genomicsDBType = isGT ? "int" : formatHeaderLine.getType().toString();
                String genomicsDBLength = isGT ? "PP" : (formatHeaderLine.getType() == VCFHeaderLineType.String)
                        ? "VAR" : getLength(formatHeaderLine);
                lengthDescriptorComponentBuilder.setVariableLengthDescriptor(genomicsDBLength);

                infoBuilder
                        .setName(formatHeaderLine.getID())
                        .addType(genomicsDBType)
                        .addLength(lengthDescriptorComponentBuilder.build());

                if (formatHeaderLine.getID().equals("DP") && dpIndex != -1) {
                    GenomicsDBVidMapProto.InfoField prevDPField = remove(infoFields, dpIndex);

                    infoBuilder
                            .addVcfFieldClass(prevDPField.getVcfFieldClass(0))
                            .addVcfFieldClass("FORMAT");
                } else {
                    infoBuilder
                            .addVcfFieldClass("FORMAT");
                }

                GenomicsDBVidMapProto.InfoField formatField = infoBuilder.build();
                infoFields.add(formatField);
                if (formatHeaderLine.getID().equals("DP")) {
                    dpIndex = infoFields.indexOf(formatField);
                }
            } else if (headerLine instanceof VCFInfoHeaderLine) {
                VCFInfoHeaderLine infoHeaderLine = (VCFInfoHeaderLine) headerLine;
                final String infoFieldName = infoHeaderLine.getID();
                infoBuilder
                        .setName(infoFieldName);
                //allele specific annotations
                if (mRLengthHistogramFieldsWithFloatBins.contains(infoFieldName)
                        || mRLengthTwoDFloatVectorFields.contains(infoFieldName)
                        || mRLengthTwoDIntVectorFields.contains(infoFieldName)
                        ) {
                    lengthDescriptorComponentBuilder.setVariableLengthDescriptor("R");
                    infoBuilder.addLength(lengthDescriptorComponentBuilder.build());
                    lengthDescriptorComponentBuilder.setVariableLengthDescriptor("var"); //ignored - can set anything here
                    infoBuilder.addLength(lengthDescriptorComponentBuilder.build());
                    infoBuilder.addVcfDelimiter("|");
                    infoBuilder.addVcfDelimiter(",");
                    if (mRLengthHistogramFieldsWithFloatBins.contains(infoFieldName)) {
                        //Each element of the vector is a tuple <float, int>
                        infoBuilder.addType("float");
                        infoBuilder.addType("int");
                        infoBuilder.setVCFFieldCombineOperation("histogram_sum");
                    } else {
                        infoBuilder.setVCFFieldCombineOperation("element_wise_sum");
                        if (mRLengthTwoDFloatVectorFields.contains(infoFieldName)) {
                            infoBuilder.addType("float");
                        } else if (mRLengthTwoDIntVectorFields.contains(infoFieldName)) {
                            infoBuilder.addType("int");
                        }
                    }
                } else {
                    lengthDescriptorComponentBuilder.setVariableLengthDescriptor(
                            infoHeaderLine.getType() == VCFHeaderLineType.String ? "var" : getLength(infoHeaderLine));

                    infoBuilder.addType(infoHeaderLine.getType().toString())
                            .addLength(lengthDescriptorComponentBuilder.build());
                }

                if (infoHeaderLine.getID().equals("DP") && dpIndex != -1) {
                    GenomicsDBVidMapProto.InfoField prevDPield = remove(infoFields, dpIndex);

                    infoBuilder
                            .addVcfFieldClass(prevDPield.getVcfFieldClass(0))
                            .addVcfFieldClass("INFO");
                } else {
                    infoBuilder
                            .addVcfFieldClass("INFO");
                }

                GenomicsDBVidMapProto.InfoField infoField = infoBuilder.build();
                infoFields.add(infoField);
                if (infoHeaderLine.getID().equals("DP")) {
                    dpIndex = infoFields.indexOf(infoField);
                }
            } else if (headerLine instanceof VCFFilterHeaderLine) {
                VCFFilterHeaderLine filterHeaderLine = (VCFFilterHeaderLine) headerLine;

                infoBuilder.setName(filterHeaderLine.getID());

                if (!filterHeaderLine.getValue().isEmpty()) {
                    infoBuilder.addType(filterHeaderLine.getValue());
                } else {
                    infoBuilder.addType("int");
                }
                infoBuilder.addVcfFieldClass("FILTER");
                GenomicsDBVidMapProto.InfoField filterField = infoBuilder.build();

                infoFields.add(filterField);
            } else if (headerLine instanceof VCFContigHeaderLine) {
                VCFContigHeaderLine contigHeaderLine = (VCFContigHeaderLine) headerLine;

                long length = contigHeaderLine.getSAMSequenceRecord().getSequenceLength();
                GenomicsDBVidMapProto.Chromosome.Builder contigBuilder =
                        GenomicsDBVidMapProto.Chromosome.newBuilder();
                contigBuilder
                        .setName(contigHeaderLine.getID())
                        .setLength(length)
                        .setTiledbColumnOffset(columnOffset);

                columnOffset += length;

                GenomicsDBVidMapProto.Chromosome chromosome = contigBuilder.build();

                contigs.add(chromosome);
            }
        }

        GenomicsDBVidMapProto.VidMappingPB.Builder vidMapBuilder =
                GenomicsDBVidMapProto.VidMappingPB.newBuilder();

        return vidMapBuilder
                .addAllFields(infoFields)
                .addAllContigs(contigs)
                .build();
    }

    private static GenomicsDBVidMapProto.InfoField remove(
            List<GenomicsDBVidMapProto.InfoField> infoFields,
            int dpIndex) {
        GenomicsDBVidMapProto.InfoField dpFormatField = infoFields.get(dpIndex);
        infoFields.remove(dpIndex);
        return dpFormatField;
    }

    /**
     * Maps the "Number" from INFO or FORMAT fields in VCF header
     * to GenomicsDB lengths. If unbounded, the field becomes a
     * variable length attribute (discussed in more detail in TileDB
     * tutorial tiledb.org). If "A", "R" or "G", the field also
     * becomes a variable length attribute, but the order and length
     * depends on order/number of alleles or genotypes. Integers
     * remain integers
     *
     * @param headerLine Info or Format header line from VCF
     * @return VAR, A, R, G, or integer values from VCF
     */
    private static String getLength(VCFHeaderLine headerLine) {

        VCFHeaderLineCount type;

        if (headerLine instanceof VCFFormatHeaderLine) {
            type = ((VCFFormatHeaderLine) headerLine).getCountType();
        } else {
            type = ((VCFInfoHeaderLine) headerLine).getCountType();
        }

        String length = "";
        int count;
        boolean isFlagType;
        switch (type) {
            case UNBOUNDED:
                length = "VAR";
                break;
            case A:
                length = "A";
                break;
            case R:
                length = "R";
                break;
            case G:
                length = "G";
                break;
            case INTEGER: {
                if (headerLine instanceof VCFFormatHeaderLine) {
                    VCFFormatHeaderLine formatHeaderLine = (VCFFormatHeaderLine) headerLine;
                    count = formatHeaderLine.getCount();
                    isFlagType = formatHeaderLine.getType().equals(VCFHeaderLineType.Flag);
                } else {
                    VCFInfoHeaderLine infoHeaderLine = (VCFInfoHeaderLine) headerLine;
                    count = infoHeaderLine.getCount();
                    isFlagType = infoHeaderLine.getType().equals(VCFHeaderLineType.Flag);
                }
                //Weird Flag fields - Number=0 in the VCF header :(
                if (count == 0 && isFlagType)
                    length = "1";
                else
                    length = String.valueOf(count);
                break;
            }
        }
        return length;
    }

    /**
     * Static function that reads sample names from the vcfHeader and adds entries to the map.
     * The function assumes that the samples will be assigned row indexes beginning at rowIdx
     * and that the sample names specified in the header
     * are globally unique (across all streams/files)
     *
     * @param sampleIndexToInfo map: key=sampleIndex in vcfHeader: value=SampleInfo
     * @param vcfHeader         VCF header
     * @param rowIdx            Starting row index from which to assign
     * @return rowIdx+#samples in the header
     */
    public static long initializeSampleInfoMapFromHeader(Map<Integer, SampleInfo> sampleIndexToInfo,
                                                         final VCFHeader vcfHeader,
                                                         final long rowIdx) {
        final List<String> headerSampleNames = vcfHeader.getGenotypeSamples();
        final int numSamplesInHeader = headerSampleNames.size();
        for (int i = 0; i < numSamplesInHeader; ++i)
            sampleIndexToInfo.put(i, new SampleInfo(headerSampleNames.get(i), rowIdx + i));
        return rowIdx + numSamplesInHeader;
    }

    /**
     * Utility function that returns a list of ChromosomeInterval objects for
     * the column partition specified by the loader JSON file and rank/partition index
     *
     * @param loaderJSONFile path to loader JSON file
     * @param partitionIdx   rank/partition index
     * @return list of ChromosomeInterval objects for the specified partition
     * @throws ParseException when there is a bug in the JNI interface and a faulty JSON is returned
     */
    private static ArrayList<ChromosomeInterval> getChromosomeIntervalsForColumnPartition(
            final String loaderJSONFile, final int partitionIdx) throws ParseException {
        final String chromosomeIntervalsJSONString =
                jniGetChromosomeIntervalsForColumnPartition(loaderJSONFile, partitionIdx);
    /* JSON format
      {
        "contigs": [
           { "chr1": [ 100, 200] },
           { "chr2": [ 500, 600] }
        ]
      }
    */
        ArrayList<ChromosomeInterval> chromosomeIntervals = new ArrayList<>();
        JSONParser parser = new JSONParser();
        JSONObject topObj = (JSONObject) (parser.parse(chromosomeIntervalsJSONString));
        assert topObj.containsKey("contigs");
        JSONArray listOfDictionaries = (JSONArray) (topObj.get("contigs"));
        for (Object currDictObj : listOfDictionaries) {
            JSONObject currDict = (JSONObject) currDictObj;
            assert currDict.size() == 1; //1 entry
            for (Object currEntryObj : currDict.entrySet()) {
                Map.Entry<String, JSONArray> currEntry = (Map.Entry<String, JSONArray>) currEntryObj;
                JSONArray currValue = currEntry.getValue();
                assert currValue.size() == 2;
                chromosomeIntervals.add(new ChromosomeInterval(currEntry.getKey(), (Long) (currValue.get(0)),
                        (Long) (currValue.get(1))));
            }
        }
        return chromosomeIntervals;
    }

    /**
     * Utility function that returns a MultiChromosomeIterator given an AbstractFeatureReader
     * that will iterate over the VariantContext objects provided by the reader belonging
     * to the column partition specified by the loader JSON file and rank/partition index
     *
     * @param <SOURCE>       LineIterator for VCFs, PositionalBufferedStream for BCFs
     * @param reader         AbstractFeatureReader over VariantContext objects -
     *                       SOURCE can vary - BCF v/s VCF for example
     * @param loaderJSONFile path to loader JSON file
     * @param partitionIdx   rank/partition index
     * @return MultiChromosomeIterator that iterates over VariantContext objects in the reader
     * belonging to the specified column partition
     * @throws IOException    when the reader's query method throws an IOException
     * @throws ParseException when there is a bug in the JNI interface and a faulty JSON is returned
     */
    public static <SOURCE> MultiChromosomeIterator<SOURCE> columnPartitionIterator(
            AbstractFeatureReader<VariantContext, SOURCE> reader,
            final String loaderJSONFile,
            final int partitionIdx) throws ParseException, IOException {
        return new MultiChromosomeIterator<SOURCE>(reader,
                GenomicsDBImporter.getChromosomeIntervalsForColumnPartition(
                        loaderJSONFile, partitionIdx));
    }

    /**
     * Creates GenomicsDBImporter object when importing VCF files (no streams)
     *
     * @param loaderJSONFile Path to loader JSON file
     * @param rank           Rank of object - corresponds to the partition index in the loader
     *                       for which this object will import data
     * @param lbRowIdx       Smallest row idx which should be imported by this object
     * @param ubRowIdx       Largest row idx which should be imported by this object
     * @return status - 0 if everything was ok, -1 otherwise
     */
    private native int jniGenomicsDBImporter(String loaderJSONFile,
                                             int rank,
                                             long lbRowIdx,
                                             long ubRowIdx);

    /**
     * Creates GenomicsDBImporter object when importing VCF files (no streams)
     *
     * @param loaderJSONFile Path to loader JSON file
     * @param rank           Rank of object - corresponds to the partition index in the
     *                       loader for which this object will import data
     * @param lbRowIdx       Smallest row idx which should be imported by this object
     * @param ubRowIdx       Largest row idx which should be imported by this object
     * @return "pointer"/address to GenomicsDBImporter object in memory,
     * if 0, then something went wrong
     */
    private native long jniInitializeGenomicsDBImporterObject(String loaderJSONFile,
                                                              int rank,
                                                              long lbRowIdx,
                                                              long ubRowIdx);

    /**
     * Copy the vid map protocol buffer to C++ through JNI
     *
     * @param genomicsDBImporterHandle Reference to a C++ GenomicsDBImporter object
     * @param vidMapAsByteArray        INFO, FORMAT, FILTER header lines and contig positions
     * @return Reference to a C++ GenomicsDBImporter object as long
     */
    private native long jniCopyVidMap(long genomicsDBImporterHandle,
                                      byte[] vidMapAsByteArray);

    /**
     * Copy the callset map protocol buffer to C++ through JNI
     *
     * @param genomicsDBImporterHandle Reference to a C++ GenomicsDBImporter object
     * @param callsetMapAsByteArray    Callset name and row index map
     * @return Reference to a C++ GenomicsDBImporter object as long
     */
    private native long jniCopyCallsetMap(long genomicsDBImporterHandle,
                                          byte[] callsetMapAsByteArray);

    /**
     * Notify importer object that a new stream is to be added
     *
     * @param genomicsDBImporterHandle "pointer" returned by jniInitializeGenomicsDBImporterObject
     * @param streamName               name of the stream
     * @param isBCF                    use BCF format to pass data to C++ layer
     * @param bufferCapacity           in bytes
     * @param buffer                   initialization buffer containing the VCF/BCF header
     * @param numValidBytesInBuffer    num valid bytes in the buffer (length of the header)
     */
    private native void jniAddBufferStream(long genomicsDBImporterHandle,
                                           String streamName,
                                           boolean isBCF,
                                           long bufferCapacity,
                                           byte[] buffer,
                                           long numValidBytesInBuffer);

    /**
     * Setup loader after all the buffer streams are added
     *
     * @param genomicsDBImporterHandle "pointer" returned by jniInitializeGenomicsDBImporterObject
     * @param callsetMappingJSON       JSON formatted string containing globally consistent callset
     *                                 name to row index mapping
     * @param usingVidMappingProtoBuf  use protocol buffer based header
     * @return maximum number of buffer stream identifiers that can be returned in
     * mExhaustedBufferStreamIdentifiers later
     * (this depends on the number of partitions and the number of buffer streams)
     */
    private native long jniSetupGenomicsDBLoader(long genomicsDBImporterHandle,
                                                 final String callsetMappingJSON,
                                                 boolean usingVidMappingProtoBuf);

    /**
     * @param handle                "pointer" returned by jniInitializeGenomicsDBImporterObject
     * @param streamIdx             stream index
     * @param partitionIdx          partition index (unused now)
     * @param buffer                buffer containing data
     * @param numValidBytesInBuffer num valid bytes in the buffer
     */
    private native void jniWriteDataToBufferStream(long handle,
                                                   int streamIdx,
                                                   int partitionIdx,
                                                   byte[] buffer,
                                                   long numValidBytesInBuffer);

    /**
     * Import the next batch of data into TileDB/GenomicsDB
     *
     * @param genomicsDBImporterHandle   "pointer" returned by jniInitializeGenomicsDBImporterObject
     * @param exhaustedBufferIdentifiers contains the list of exhausted buffer stream identifiers
     *                                   - the number of
     *                                   exhausted streams is stored in the last element of the array
     * @return true if the whole import process is completed, false otherwise
     */
    private native boolean jniImportBatch(long genomicsDBImporterHandle,
                                          long[] exhaustedBufferIdentifiers);

    private GenomicsDBImportConfiguration.ImportConfiguration createImportConfiguration(
            String workspace,
            String arrayname,
            Long sizePerColumnPartition,
            Long segmentSize,
            boolean failIfUpdating) {

        String name = (arrayname.isEmpty()) ? DEFAULT_ARRAY_NAME : arrayname;

        GenomicsDBImportConfiguration.Partition.Builder pB = GenomicsDBImportConfiguration.Partition.newBuilder();
        GenomicsDBImportConfiguration.Partition p0 = pB.setBegin(0).setWorkspace(workspace).setArray(name).build();
        GenomicsDBImportConfiguration.ImportConfiguration.Builder importBuilder =
                GenomicsDBImportConfiguration.ImportConfiguration.newBuilder();

        return importBuilder
                .setRowBasedPartitioning(false)
                .setSizePerColumnPartition(sizePerColumnPartition)
                .addColumnPartitions(p0)
                .setProduceTiledbArray(true)
                .setNumCellsPerTile(DEFAULT_TILEDB_CELLS_PER_TILE)
                .setCompressTiledbArray(true)
                .setSegmentSize(segmentSize)
                .setTreatDeletionsAsIntervals(true)
                .setFailIfUpdating(failIfUpdating)
                .build();
    }

    /**
     * Initialize variables
     *
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank           Rank of this process (TileDB/GenomicsDB partition idx)
     * @param lbRowIdx       Smallest row idx which should be imported by this object
     * @param ubRowIdx       Largest row idx which should be imported by this object
     */
    private void initialize(String loaderJSONFile, int rank, long lbRowIdx, long ubRowIdx) {
        mLoaderJSONFile = loaderJSONFile;
        mRank = rank;
        mLbRowIdx = lbRowIdx;
        mUbRowIdx = ubRowIdx;
    }

    /**
     * Add a buffer stream as the data source - caller must:
     * 1. Call setupGenomicsDBImporter() once all streams are added
     * 2. Provide VC objects using the add() function
     * 3. Call importBatch()
     * 4. Get list of exhausted buffer streams using getNumExhaustedBufferStreams(),
     * getExhaustedBufferStreamIndex()
     * 5. If !isDone() goto 2
     *
     * @param streamName        Name of the stream being added - must be unique with respect to this
     *                          GenomicsDBImporter object
     * @param vcfHeader         VCF header for the stream
     * @param bufferCapacity    Capacity of the stream buffer in bytes
     * @param streamType        BCF_STREAM or VCF_STREAM
     * @param sampleIndexToInfo map from sample index in the vcfHeader to SampleInfo object which
     *                          contains row index and globally unique name
     *                          can be set to null, which implies that the mapping is stored in a callsets JSON file
     * @return returns buffer stream index
     */
    public int addBufferStream(final String streamName,
                               final VCFHeader vcfHeader,
                               final long bufferCapacity,
                               final VariantContextWriterBuilder.OutputType streamType,
                               final Map<Integer, SampleInfo> sampleIndexToInfo)
            throws GenomicsDBException {
        return addBufferStream(streamName, vcfHeader, bufferCapacity,
                streamType, null, sampleIndexToInfo);
    }

    /**
     * Add a sorted VC iterator as the data source - caller must:
     * 1. Call setupGenomicsDBImporter() once all iterators are added
     * 2. Call importBatch()
     * 3. Done!
     *
     * @param streamName        Name of the stream being added - must be unique with respect
     *                          to this GenomicsDBImporter object
     * @param vcfHeader         VCF header for the stream
     * @param vcIterator        Iterator over VariantContext objects
     * @param bufferCapacity    Capacity of the stream buffer in bytes
     * @param streamType        BCF_STREAM or VCF_STREAM
     * @param sampleIndexToInfo map from sample index in the vcfHeader to SampleInfo object
     *                          which contains row index and globally unique name
     *                          can be set to null, which implies that the mapping is
     *                          stored in a callsets JSON file
     * @return returns the stream index
     */
    public int addSortedVariantContextIterator(
            final String streamName,
            final VCFHeader vcfHeader,
            Iterator<VariantContext> vcIterator,
            final long bufferCapacity,
            final VariantContextWriterBuilder.OutputType streamType,
            final Map<Integer, SampleInfo> sampleIndexToInfo) throws GenomicsDBException {
        return addBufferStream(streamName, vcfHeader, bufferCapacity, streamType,
                vcIterator, sampleIndexToInfo);
    }

    /**
     * Sets sorted VC iterator as the data source and calls setupGenomicsDBImporter().
     * No more streams/iterators can
     * be added after this function is called. The caller must:
     * 1. Call importBatch()
     * 2. Done!
     *
     * @param streamName        Name of the stream being added - must be unique with respect to
     *                          this GenomicsDBImporter object
     * @param vcfHeader         VCF header for the stream
     * @param vcIterator        Iterator over VariantContext objects
     * @param bufferCapacity    Capacity of the stream buffer in bytes
     * @param streamType        BCF_STREAM or VCF_STREAM
     * @param sampleIndexToInfo map from sample index in the vcfHeader to SampleInfo
     *                          object which contains row index and globally unique name
     *                          can be set to null, which implies that the mapping is stored in a
     *                          callsets JSON file
     * @return returns the stream index
     * @throws GenomicsDBException thrown if incorrect iterator or missing JSON configuration
     * @throws IOException         thrown if incorrect iterator or missing JSON configuration
     *                             files
     */
    public int setSortedVariantContextIterator(
            final String streamName,
            final VCFHeader vcfHeader,
            Iterator<VariantContext> vcIterator,
            final long bufferCapacity,
            final VariantContextWriterBuilder.OutputType streamType,
            final Map<Integer, SampleInfo> sampleIndexToInfo) throws GenomicsDBException, IOException {
        int streamIdx = addSortedVariantContextIterator(streamName, vcfHeader, vcIterator,
                bufferCapacity, streamType, sampleIndexToInfo);
        setupGenomicsDBImporter();
        return streamIdx;
    }

    /**
     * Add a buffer stream or VC iterator - internal function
     *
     * @param streamName        Name of the stream being added - must be unique with respect to this
     *                          GenomicsDBImporter object
     * @param vcfHeader         VCF header for the stream
     * @param bufferCapacity    Capacity of the stream buffer in bytes
     * @param streamType        BCF_STREAM or VCF_STREAM
     * @param vcIterator        Iterator over VariantContext objects - can be null
     * @param sampleIndexToInfo map from sample index in the vcfHeader to SampleInfo object which
     *                          contains row index and globally unique name can be set to null,
     *                          which implies that the mapping is stored in a callsets JSON file
     * @return returns the stream index
     */
    @SuppressWarnings("unchecked")
    private int addBufferStream(final String streamName,
                                final VCFHeader vcfHeader,
                                final long bufferCapacity,
                                final VariantContextWriterBuilder.OutputType streamType,
                                Iterator<VariantContext> vcIterator,
                                final Map<Integer, SampleInfo> sampleIndexToInfo)
            throws GenomicsDBException {
        if (mIsLoaderSetupDone)
            throw new GenomicsDBException("Cannot add buffer streams after "
                    + "setupGenomicsDBImporter() is called");
        //First time a buffer is added
        if (!mContainsBufferStreams) {
            if (mGenomicsDBImporterObjectHandle == 0) {
                mGenomicsDBImporterObjectHandle = jniInitializeGenomicsDBImporterObject(mLoaderJSONFile,
                        mRank, mLbRowIdx, mUbRowIdx);
            }
            if (mGenomicsDBImporterObjectHandle == 0)
                throw new GenomicsDBException("Could not initialize GenomicsDBImporter object");
            mBufferStreamWrapperVector = new ArrayList<>();
            mCallsetMappingJSON = new JSONObject();
            mContainsBufferStreams = true;
        }
        mBufferStreamWrapperVector.add(new GenomicsDBImporterStreamWrapper(vcfHeader,
                bufferCapacity, streamType, vcIterator));
        int currIdx = mBufferStreamWrapperVector.size() - 1;
        SilentByteBufferStream currStream = mBufferStreamWrapperVector.get(currIdx).mStream;
        jniAddBufferStream(mGenomicsDBImporterObjectHandle,
                streamName, streamType == VariantContextWriterBuilder.OutputType.BCF_STREAM,
                bufferCapacity, currStream.getBuffer(), currStream.getNumValidBytes());
        if (sampleIndexToInfo != null) {
            for (Map.Entry<Integer, SampleInfo> currEntry : sampleIndexToInfo.entrySet()) {
                JSONObject sampleJSON = new JSONObject();
                sampleJSON.put("row_idx", currEntry.getValue().mRowIdx);
                sampleJSON.put("stream_name", streamName);
                sampleJSON.put("idx_in_file", currEntry.getKey());
                mCallsetMappingJSON.put(currEntry.getValue().mName, sampleJSON);
            }
        }
        return currIdx;
    }

    /**
     * Setup the importer after all the buffer streams are added, but before any
     * data is inserted into any stream
     * No more buffer streams can be added once setupGenomicsDBImporter() is called
     *
     * @throws IOException throws IOException if modified callsets JSON cannot be written
     */
    @SuppressWarnings("unchecked")
    public void setupGenomicsDBImporter() throws IOException {
        if (mIsLoaderSetupDone)
            return;
        //Callset mapping JSON - convert to string
        JSONObject topCallsetJSON = new JSONObject();
        topCallsetJSON.put("callsets", mCallsetMappingJSON);
        StringWriter stringWriter = new StringWriter();
        topCallsetJSON.writeJSONString(stringWriter);
        //Call native setupGenomicsDBImporter()
        long mMaxBufferStreamIdentifiers = jniSetupGenomicsDBLoader(mGenomicsDBImporterObjectHandle,
                stringWriter.toString(), mUsingVidMappingProtoBuf);
        //Why 2* - each identifier is a pair<buffer_stream_idx, partition_idx>
        //Why +1 - the last element will contain the number of exhausted stream identifiers
        //when importBatch() is called
        mExhaustedBufferStreamIdentifiers = new long[2 * ((int) mMaxBufferStreamIdentifiers) + 1];
        //Set all streams to empty
        //Add all streams to mExhaustedBufferStreamIdentifiers - this way when importBatch is
        //called the first time all streams' data are written
        for (int i = 0, idx = 0; i < mBufferStreamWrapperVector.size(); ++i, idx += 2) {
            SilentByteBufferStream currStream = mBufferStreamWrapperVector.get(i).mStream;
            currStream.setNumValidBytes(0);
            mExhaustedBufferStreamIdentifiers[idx] = i;
            mExhaustedBufferStreamIdentifiers[idx + 1] = 0;
        }
        //Set number of exhausted buffer streams - all streams are exhausted the first time
        mNumExhaustedBufferStreams = mBufferStreamWrapperVector.size();
        mExhaustedBufferStreamIdentifiers[(int) (2 * mMaxBufferStreamIdentifiers)] =
                mNumExhaustedBufferStreams;
        mIsLoaderSetupDone = true;
    }

    /**
     * Write VariantContext object to stream - may fail if the buffer is full
     * It's the caller's responsibility keep track of the VC object that's not written
     *
     * @param vc        VariantContext object
     * @param streamIdx index of the stream returned by the addBufferStream() call
     * @return true if the vc object was written successfully, false otherwise
     */
    public boolean add(VariantContext vc, final int streamIdx)
            throws GenomicsDBException, RuntimeIOException {
        if (streamIdx < 0 || streamIdx >= mBufferStreamWrapperVector.size())
            throw new GenomicsDBException("Invalid stream idx "
                    + Integer.toString(streamIdx) + " must be between [0-"
                    + Long.toString(mBufferStreamWrapperVector.size() - 1) + "]");
        if (!mIsLoaderSetupDone)
            throw new GenomicsDBException("Cannot add VariantContext objects " +
                    "to streams before calling setupGenomicsDBImporter()");
        GenomicsDBImporterStreamWrapper currWrapper = mBufferStreamWrapperVector.get(streamIdx);
        currWrapper.mVCWriter.add(vc);
        SilentByteBufferStream currStream = currWrapper.mStream;
        //at least one record already existed in the buffer
        if (currStream.overflow() && currStream.getMarker() > 0) {
            //Set num valid bytes to marker - marker points to location
            //after the last valid serialized vc in the buffer
            currStream.setNumValidBytes(currStream.getMarker());
            return false;
        } else {
            //the first record to be added to the buffer is too large, resize buffer
            while (currStream.overflow()) {
                currStream.resize(2 * currStream.size() + 1);
                currStream.setNumValidBytes(0);
                currStream.setOverflow(false);
                currWrapper.mVCWriter.add(vc);
            }
            //Update marker - marker points to location after the last valid serialized vc in the buffer
            currStream.setMarker(currStream.getNumValidBytes());
            return true;
        }
    }

    /**
     * @return true if the import process is done
     * @throws IOException if the wimport fails
     */
    public boolean importBatch() throws IOException {
        if (mDone)
            return true;
        if (!mIsLoaderSetupDone)
            setupGenomicsDBImporter();
        boolean allExhaustedStreamsHaveIterators = true;
        while (!mDone && allExhaustedStreamsHaveIterators) {
            //Write data from buffer streams exhausted in the previous round into GenomicsDB
            for (int i = 0, idx = 0; i < mNumExhaustedBufferStreams; ++i, idx += 2) {
                int bufferStreamIdx = (int) mExhaustedBufferStreamIdentifiers[idx];
                GenomicsDBImporterStreamWrapper currWrapper =
                        mBufferStreamWrapperVector.get(bufferStreamIdx);
                //If iterator is provided, get data from iterator
                if (currWrapper.hasIterator()) {
                    while (currWrapper.getCurrentVC() != null) {
                        boolean added = add(currWrapper.getCurrentVC(), bufferStreamIdx);
                        if (added)
                            currWrapper.next();
                        else
                            break; //buffer full
                    }
                }
                SilentByteBufferStream currStream = currWrapper.mStream;
                jniWriteDataToBufferStream(mGenomicsDBImporterObjectHandle, bufferStreamIdx,
                        0, currStream.getBuffer(), currStream.getNumValidBytes());
            }
            mDone = jniImportBatch(mGenomicsDBImporterObjectHandle, mExhaustedBufferStreamIdentifiers);
            mNumExhaustedBufferStreams =
                    mExhaustedBufferStreamIdentifiers[mExhaustedBufferStreamIdentifiers.length - 1];
            //Reset markers, numValidBytesInBuffer and overflow flag for the exhausted streams
            for (long i = 0, idx = 0; i < mNumExhaustedBufferStreams; ++i, idx += 2) {
                int bufferStreamIdx = (int) mExhaustedBufferStreamIdentifiers[(int) idx];
                GenomicsDBImporterStreamWrapper currWrapper =
                        mBufferStreamWrapperVector.get(bufferStreamIdx);
                if (!currWrapper.hasIterator())
                    allExhaustedStreamsHaveIterators = false;
                SilentByteBufferStream currStream = currWrapper.mStream;
                currStream.setOverflow(false);
                currStream.setMarker(0);
                currStream.setNumValidBytes(0);
            }
            if (mDone) {
                mGenomicsDBImporterObjectHandle = 0;
                mContainsBufferStreams = false;
                mIsLoaderSetupDone = false;
            }
        }

        return mDone;
    }

    /**
     * @return get number of buffer streams for which new data must be supplied
     */
    public long getNumExhaustedBufferStreams() {
        return mNumExhaustedBufferStreams;
    }

    /**
     * Get buffer stream index of i-th exhausted stream
     * There are mNumExhaustedBufferStreams and the caller must provide data for streams
     * with indexes getExhaustedBufferStreamIndex(0), getExhaustedBufferStreamIndex(1),...,
     * getExhaustedBufferStreamIndex(mNumExhaustedBufferStreams-1)
     *
     * @param i i-th exhausted buffer stream
     * @return the buffer stream index of the i-th exhausted stream
     */
    public int getExhaustedBufferStreamIndex(final long i) {
        assert i < mNumExhaustedBufferStreams && i >= 0;
        //why 2* - exhausted buffer stream identifier is a pair<stream_idx, partition_idx>
        return (int) mExhaustedBufferStreamIdentifiers[2 * ((int) i)];
    }

    /**
     * Is the import process completed
     *
     * @return true if complete, false otherwise
     */
    public boolean isDone() {
        return mDone;
    }

    /**
     * Utility function that returns a MultiChromosomeIterator given an AbstractFeatureReader
     * that will iterate over the VariantContext objects provided by the reader belonging
     * to the column partition specified by this object's loader JSON file and rank/partition index
     *
     * @param <SOURCE> LineIterator for VCFs, PositionalBufferedStream for BCFs
     * @param reader   AbstractFeatureReader over VariantContext objects -
     *                 SOURCE can vary - BCF v/s VCF for example
     * @return MultiChromosomeIterator that iterates over VariantContext objects in the reader
     * belonging to the specified column partition
     * @throws IOException    when the reader's query method throws an IOException
     * @throws ParseException when there is a bug in the JNI interface and a faulty JSON is returned
     */
    public <SOURCE> MultiChromosomeIterator<SOURCE> columnPartitionIterator(
            AbstractFeatureReader<VariantContext, SOURCE> reader) throws ParseException, IOException {
        return GenomicsDBImporter.columnPartitionIterator(reader, mLoaderJSONFile, mRank);
    }

    /**
     * Write to TileDB/GenomicsDB using the configuration specified in the
     * loader file passed to constructor
     */
    public void write() throws GenomicsDBException {
        write(mLoaderJSONFile, mRank, 0, Long.MAX_VALUE - 1);
    }

    /**
     * Write to TileDB/GenomicsDB using the configuration specified in the
     * loader file passed to constructor
     *
     * @param lbRowIdx Minimum row idx from which new data will be added
     */
    public void write(long lbRowIdx) throws GenomicsDBException {
        write(mLoaderJSONFile, mRank, lbRowIdx, Long.MAX_VALUE - 1);
    }

    /**
     * Write to TileDB/GenomicsDB using the configuration specified in the
     * loader file passed to constructor
     *
     * @param rank     Rank of this process (TileDB/GenomicsDB partition idx)
     * @param lbRowIdx Minimum row idx from which new data will be added
     */
    public void write(int rank, long lbRowIdx) throws GenomicsDBException {
        write(mLoaderJSONFile, rank, lbRowIdx, Long.MAX_VALUE - 1);
    }

    /**
     * Write to TileDB/GenomicsDB using the configuration specified in the
     * loader file passed to constructor
     *
     * @param rank     Rank of this process (TileDB/GenomicsDB partition idx)
     * @param lbRowIdx Minimum row idx from which new data will be added
     * @param ubRowIdx Maximum row idx upto which new data will be added
     */
    public void write(int rank, long lbRowIdx, long ubRowIdx) throws GenomicsDBException {
        write(mLoaderJSONFile, rank, lbRowIdx, ubRowIdx);
    }

    /**
     * Write to TileDB/GenomicsDB
     *
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank           Rank of this process (TileDB/GenomicsDB partition idx)
     * @param lbRowIdx       Minimum row idx from which new data will be added
     * @param ubRowIdx       Maximum row idx upto which new data will be added
     */
    public void write(String loaderJSONFile, int rank, long lbRowIdx, long ubRowIdx)
            throws GenomicsDBException {
        mDone = false;
        if (loaderJSONFile == null)
            throw new GenomicsDBException("Loader JSON file not specified");
        if (mContainsBufferStreams)
            throw new GenomicsDBException("Cannot call write() functions if buffer streams are added");
        int status = jniGenomicsDBImporter(loaderJSONFile, rank, lbRowIdx, ubRowIdx);
        if (status != 0)
            throw new GenomicsDBException("GenomicsDBImporter write failed for loader JSON: "
                    + loaderJSONFile + " rank: " + rank);
        mDone = true;
    }
}
