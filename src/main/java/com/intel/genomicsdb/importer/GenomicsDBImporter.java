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

import com.intel.genomicsdb.GenomicsDBLibLoader;
import com.intel.genomicsdb.exception.GenomicsDBException;
import com.intel.genomicsdb.importer.extensions.CallSetMapExtensions;
import com.intel.genomicsdb.importer.extensions.JsonFileExtensions;
import com.intel.genomicsdb.importer.extensions.VidMapExtensions;
import com.intel.genomicsdb.importer.model.ChromosomeInterval;
import com.intel.genomicsdb.importer.model.SampleInfo;
import com.intel.genomicsdb.model.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static com.intel.genomicsdb.Constants.CHROMOSOME_INTERVAL_FOLDER;
import static com.intel.genomicsdb.importer.Constants.*;

/**
 * Java wrapper for vcf2tiledb - imports VCFs into TileDB/GenomicsDB.
 * All vid information is assumed to be set correctly by the user (JSON files)
 */
public class GenomicsDBImporter extends GenomicsDBImporterJni implements JsonFileExtensions, CallSetMapExtensions,
        VidMapExtensions {
    static {
        try {
            boolean loaded = GenomicsDBLibLoader.loadLibrary();
            if (!loaded) throw new GenomicsDBException("Could not load genomicsdb native library");
        } catch (UnsatisfiedLinkError ule) {
            throw new GenomicsDBException("Could not load genomicsdb native library", ule);
        }
    }

    private ImportConfig config;
    private GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB;
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
     * Constructor
     *
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     */
    public GenomicsDBImporter(final String loaderJSONFile) {
        initialize(loaderJSONFile, 0, 0, Long.MAX_VALUE - 1);
    }

    /**
     * Constructor
     *
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank           Rank of this process (TileDB/GenomicsDB partition idx)
     */
    public GenomicsDBImporter(final String loaderJSONFile, final int rank) {
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
    public GenomicsDBImporter(final String loaderJSONFile, final int rank, final long lbRowIdx, final long ubRowIdx) {
        initialize(loaderJSONFile, rank, lbRowIdx, ubRowIdx);
    }

    /**
     * Constructor to create required data structures from a list
     * of GVCF files and a chromosome interval. This constructor
     * is developed specifically for GATK4 GenomicsDBImport tool.
     *
     * @param sampleToReaderMap Feature Readers objects corresponding to input GVCF files
     * @param rank              Rank of object - corresponds to the partition index in the loader
     * @throws IOException when load into TileDB array fails
     */
    private GenomicsDBImporter(final ImportConfig importConfig,
                               final Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
                               final int rank) throws IOException {
        // Mark this flag so that protocol buffer based vid
        // and callset map are propagated to C++ GenomicsDBImporter
        mUsingVidMappingProtoBuf = true;

        GenomicsDBVidMapProto.VidMappingPB vidMapPB = generateVidMapFromMergedHeader(importConfig.getMergedHeader());

        callsetMappingPB = generateSortedCallSetMap(sampleToReaderMap, importConfig.isUseSamplesInOrder(),
                importConfig.isValidateSampleToReaderMap(), importConfig.getImportConfiguration()
                        .getLbCallsetRowIdx());

        File importJSONFile = dumpTemporaryLoaderJSONFile(importConfig.getImportConfiguration(), "");

        initialize(importJSONFile.getAbsolutePath(), rank,
                importConfig.getImportConfiguration().getLbCallsetRowIdx(),
                importConfig.getImportConfiguration().getUbCallsetRowIdx());

        mGenomicsDBImporterObjectHandle = jniInitializeGenomicsDBImporterObject(mLoaderJSONFile, mRank,
                importConfig.getImportConfiguration().getLbCallsetRowIdx(),
                importConfig.getImportConfiguration().getUbCallsetRowIdx());

        jniCopyVidMap(mGenomicsDBImporterObjectHandle, vidMapPB.toByteArray());
        jniCopyCallsetMap(mGenomicsDBImporterObjectHandle, callsetMappingPB.toByteArray());

        String chromosomeName = importConfig.getImportConfiguration().getColumnPartitions(rank).getBegin().getContigPosition().getContig();
        int chromosomeStart = (int) importConfig.getImportConfiguration().getColumnPartitions(rank).getBegin().getContigPosition().getPosition();
        int chromosomeEnd = (int) importConfig.getImportConfiguration().getColumnPartitions(rank).getEnd().getContigPosition().getPosition();

        for (GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap sampleToIDMap :
                callsetMappingPB.getCallsetsList()) {

            String sampleName = sampleToIDMap.getSampleName();

            FeatureReader<VariantContext> featureReader = sampleToReaderMap.get(sampleName);

            CloseableIterator<VariantContext> iterator = featureReader.query(chromosomeName, chromosomeStart, chromosomeEnd);

            addSortedVariantContextIterator(
                    sampleToIDMap.getStreamName(),
                    (VCFHeader) featureReader.getHeader(),
                    iterator,
                    importConfig.getImportConfiguration().getSizePerColumnPartition(),
                    importConfig.isPassAsVcf() ? VariantContextWriterBuilder.OutputType.VCF_STREAM
                            : VariantContextWriterBuilder.OutputType.BCF_STREAM,
                    null);
        }
    }

    /**
     * Constructor to create required data structures from a list
     * of GVCF files and a chromosome interval. This constructor
     * is developed specifically for running Chromosome intervals imports in parallel.
     *
     * @param config Parallel import configuration
     * @throws FileNotFoundException
     */
    public GenomicsDBImporter(final ImportConfig config) throws FileNotFoundException {
        this.config = config;
        this.config.setImportConfiguration(addExplicitValuesToImportConfiguration(config));
        //This sorts the list sampleNames if !useSamplesInOrder
        //Why you should use this? If you are writing multiple partitions in different machines,
        //you must have consistent ordering of samples across partitions. If file order is different
        //in different processes, then set useSamplesInOrder to false and let the sort in
        //generateSortedCallSetMap ensure consistent ordering across samples
        callsetMappingPB = this.generateSortedCallSetMap(new ArrayList<>(this.config.getSampleNameToVcfPath().keySet()),
                this.config.isUseSamplesInOrder());

        //Write out callset map if needed
        String outputCallsetmapJsonFilePath = this.config.getOutputCallsetmapJsonFile();
        if (outputCallsetmapJsonFilePath != null && !outputCallsetmapJsonFilePath.isEmpty())
            this.writeCallsetMapJSONFile(outputCallsetmapJsonFilePath, callsetMappingPB);

        //Write out vidmap if needed
        String vidmapOutputFilepath = this.config.getOutputVidmapJsonFile();
        if (vidmapOutputFilepath != null && !vidmapOutputFilepath.isEmpty())
            this.writeVidMapJSONFile(vidmapOutputFilepath, generateVidMapFromMergedHeader(this.config.getMergedHeader()));

        //Write out merged header if needed
        String vcfHeaderOutputFilepath = this.config.getImportConfiguration().getColumnPartitions(0).getVcfOutputFilename();
        if (vcfHeaderOutputFilepath != null && !vcfHeaderOutputFilepath.isEmpty())
            this.writeVcfHeaderFile(vcfHeaderOutputFilepath, this.config.getMergedHeader());

        //Create workspace folder to avoid issues with concurrency
        String workspace = this.config.getImportConfiguration().getColumnPartitions(0).getWorkspace();
        if (!new File(workspace).exists()) {
            int tileDBWorkspace = createTileDBWorkspace(workspace);
            if (tileDBWorkspace < 0)
                throw new IllegalStateException(String.format("Cannot create '%s' workspace.", workspace));
        }
    }

    private GenomicsDBImportConfiguration.ImportConfiguration addExplicitValuesToImportConfiguration(ImportConfig config) {
        GenomicsDBImportConfiguration.ImportConfiguration.Builder importConfigurationBuilder =
                config.getImportConfiguration().toBuilder();
        importConfigurationBuilder.setSegmentSize(config.getImportConfiguration().getSegmentSize())
                .setFailIfUpdating(config.getImportConfiguration().getFailIfUpdating())
                //TODO: making the following attributes explicit since the C++ layer is not working with the
                // protobuf object and it's defaults
                .setTreatDeletionsAsIntervals(true).setCompressTiledbArray(true).setNumCellsPerTile(1000)
                .setRowBasedPartitioning(false).setProduceTiledbArray(true).build();
        return importConfigurationBuilder.build();
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
        final String chromosomeIntervalsJSONString = jniGetChromosomeIntervalsForColumnPartition(loaderJSONFile, partitionIdx);
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
        return new MultiChromosomeIterator<>(reader, GenomicsDBImporter.getChromosomeIntervalsForColumnPartition(
                loaderJSONFile, partitionIdx));
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
     * Add a sorted VC iterator as the data source - caller must:
     * 1. Call setupGenomicsDBImporter() once all iterators are added
     * 2. Call doSingleImport()
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
        return addBufferStream(streamName, vcfHeader, bufferCapacity, streamType, vcIterator, sampleIndexToInfo);
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
    public int addBufferStream(final String streamName,
                               final VCFHeader vcfHeader,
                               final long bufferCapacity,
                               final VariantContextWriterBuilder.OutputType streamType,
                               Iterator<VariantContext> vcIterator,
                               final Map<Integer, SampleInfo> sampleIndexToInfo)
            throws GenomicsDBException {
        if (mIsLoaderSetupDone) throw new GenomicsDBException(
                "Cannot add buffer streams after setupGenomicsDBImporter() is called");
        //First time a buffer is added
        if (!mContainsBufferStreams) {
            if (mGenomicsDBImporterObjectHandle == 0) {
                mGenomicsDBImporterObjectHandle = jniInitializeGenomicsDBImporterObject(mLoaderJSONFile, mRank,
                        mLbRowIdx, mUbRowIdx);
            }
            if (mGenomicsDBImporterObjectHandle == 0) throw new GenomicsDBException(
                    "Could not initialize GenomicsDBImporter object");
            mBufferStreamWrapperVector = new ArrayList<>();
            mCallsetMappingJSON = new JSONObject();
            mContainsBufferStreams = true;
        }
        mBufferStreamWrapperVector.add(new GenomicsDBImporterStreamWrapper(vcfHeader, bufferCapacity, streamType, vcIterator));
        int currIdx = mBufferStreamWrapperVector.size() - 1;
        SilentByteBufferStream currStream = mBufferStreamWrapperVector.get(currIdx).getStream();
        jniAddBufferStream(mGenomicsDBImporterObjectHandle, streamName, streamType == VariantContextWriterBuilder.OutputType.BCF_STREAM,
                bufferCapacity, currStream.getBuffer(), currStream.getNumValidBytes());
        if (sampleIndexToInfo != null) {
            for (Map.Entry<Integer, SampleInfo> currEntry : sampleIndexToInfo.entrySet()) {
                JSONObject sampleJSON = new JSONObject();
                sampleJSON.put("row_idx", currEntry.getValue().getRowIdx());
                sampleJSON.put("stream_name", streamName);
                sampleJSON.put("idx_in_file", currEntry.getKey());
                mCallsetMappingJSON.put(currEntry.getValue().getName(), sampleJSON);
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
        if (mIsLoaderSetupDone) return;
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
        //when doSingleImport() is called
        mExhaustedBufferStreamIdentifiers = new long[2 * ((int) mMaxBufferStreamIdentifiers) + 1];
        //Set all streams to empty
        //Add all streams to mExhaustedBufferStreamIdentifiers - this way when doSingleImport is
        //called the first time all streams' data are written
        for (int i = 0, idx = 0; i < mBufferStreamWrapperVector.size(); ++i, idx += 2) {
            SilentByteBufferStream currStream = mBufferStreamWrapperVector.get(i).getStream();
            currStream.setNumValidBytes(0);
            mExhaustedBufferStreamIdentifiers[idx] = i;
            mExhaustedBufferStreamIdentifiers[idx + 1] = 0;
        }
        //Set number of exhausted buffer streams - all streams are exhausted the first time
        mNumExhaustedBufferStreams = mBufferStreamWrapperVector.size();
        mExhaustedBufferStreamIdentifiers[(int) (2 * mMaxBufferStreamIdentifiers)] = mNumExhaustedBufferStreams;
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
    public boolean add(VariantContext vc, final int streamIdx) throws GenomicsDBException, RuntimeIOException {
        if (streamIdx < 0 || streamIdx >= mBufferStreamWrapperVector.size()) throw new GenomicsDBException(
                "Invalid stream idx " + Integer.toString(streamIdx) + " must be between [0-"
                        + Long.toString(mBufferStreamWrapperVector.size() - 1) + "]");
        if (!mIsLoaderSetupDone) throw new GenomicsDBException("Cannot add VariantContext objects to streams before " +
                "calling setupGenomicsDBImporter()");
        GenomicsDBImporterStreamWrapper currWrapper = mBufferStreamWrapperVector.get(streamIdx);
        currWrapper.getVcWriter().add(vc);
        SilentByteBufferStream currStream = currWrapper.getStream();
        //at least one record already existed in the buffer
        if (currStream.overflow() && currStream.getMarker() > 0) {
            //Set num valid bytes to marker - marker points to location
            //after the last valid serialized vc in the buffer
            currStream.setNumValidBytes(currStream.getMarker());
            return false;
        }
        //the first record to be added to the buffer is too large, resize buffer
        while (currStream.overflow()) {
            currStream.resize(2 * currStream.size() + 1);
            currStream.setNumValidBytes(0);
            currStream.setOverflow(false);
            currWrapper.getVcWriter().add(vc);
        }
        //Update marker - marker points to location after the last valid serialized vc in the buffer
        currStream.setMarker(currStream.getNumValidBytes());
        return true;
    }

    /**
     * Only to be used in cases where iterator of VariantContext are not used. The data is written to buffers
     * directly after which this function is called. See TestBufferStreamGenomicsDBImporter.java for an example
     *
     * @return true if the import process is done
     * @throws IOException if the import fails
     */
    public boolean doSingleImport() throws IOException {
        if (mDone) return true;
        if (!mIsLoaderSetupDone) setupGenomicsDBImporter();
        boolean allExhaustedStreamsHaveIterators = true;
        while (!mDone && allExhaustedStreamsHaveIterators) {
            //Write data from buffer streams exhausted in the previous round into GenomicsDB
            for (int i = 0, idx = 0; i < mNumExhaustedBufferStreams; ++i, idx += 2) {
                int bufferStreamIdx = (int) mExhaustedBufferStreamIdentifiers[idx];
                GenomicsDBImporterStreamWrapper currWrapper = mBufferStreamWrapperVector.get(bufferStreamIdx);
                //If iterator is provided, get data from iterator
                if (currWrapper.hasIterator()) {
                    while (currWrapper.getCurrentVC() != null) {
                        boolean added = add(currWrapper.getCurrentVC(), bufferStreamIdx);
                        if (added) currWrapper.next();
                        else break; //buffer full
                    }
                }
                SilentByteBufferStream currStream = currWrapper.getStream();
                jniWriteDataToBufferStream(mGenomicsDBImporterObjectHandle, bufferStreamIdx, 0, currStream.getBuffer(),
                        currStream.getNumValidBytes());
            }
            mDone = jniImportBatch(mGenomicsDBImporterObjectHandle, mExhaustedBufferStreamIdentifiers);
            mNumExhaustedBufferStreams = mExhaustedBufferStreamIdentifiers[mExhaustedBufferStreamIdentifiers.length - 1];
            //Reset markers, numValidBytesInBuffer and overflow flag for the exhausted streams
            for (long i = 0, idx = 0; i < mNumExhaustedBufferStreams; ++i, idx += 2) {
                int bufferStreamIdx = (int) mExhaustedBufferStreamIdentifiers[(int) idx];
                GenomicsDBImporterStreamWrapper currWrapper = mBufferStreamWrapperVector.get(bufferStreamIdx);
                if (!currWrapper.hasIterator()) allExhaustedStreamsHaveIterators = false;
                SilentByteBufferStream currStream = currWrapper.getStream();
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
     * Import multiple chromosome interval
     *
     * @throws InterruptedException when there is an exception in any of the threads in the stream
     */
    public void executeImport() throws InterruptedException {
        executeImport(0);
    }

    /**
     * Import multiple chromosome interval
     *
     */
    public void executeImport(final int numThreads) {
        final int batchSize = this.config.getBatchSize();
        final int sampleCount = this.config.getSampleNameToVcfPath().size();
        final int updatedBatchSize = (batchSize == DEFAULT_ZERO_BATCH_SIZE) ? sampleCount : batchSize;
        final int numberPartitions = this.config.getImportConfiguration().getColumnPartitionsList().size();

        ExecutorService executor = numThreads == 0 ? ForkJoinPool.commonPool() : Executors.newFixedThreadPool(numThreads);

        //Set size_per_column_partition once
        this.config.setImportConfiguration(this.config.getImportConfiguration().toBuilder()
            .setSizePerColumnPartition(this.config.getImportConfiguration().getSizePerColumnPartition()
              * sampleCount).build());

        //Iterate over sorted sample list in batches
        for (int i = 0, batchCount = 1; i < sampleCount; i += updatedBatchSize, ++batchCount) {
            final int index = i;

            IntStream.range(0, numberPartitions).forEach(rank ->
                    updateConfigPartitionsAndLbUb(this.config, index, rank));

            List<CompletableFuture<Boolean>> futures = IntStream.range(0, numberPartitions).mapToObj(rank ->
                    CompletableFuture.supplyAsync(() -> {
                        try {
                            final Map<String, FeatureReader<VariantContext>> sampleToReaderMap =
                                    this.config.sampleToReaderMapCreator().apply(
                                            this.config.getSampleNameToVcfPath(), updatedBatchSize, index);
                            GenomicsDBImporter importer = new GenomicsDBImporter(this.config, sampleToReaderMap, rank);
                            return importer.doSingleImport();
                        } catch (IOException ex) {
                            throw new IllegalStateException("There was an unhandled exception during chromosome interval import.", ex);
                        }
                    }, executor)
            ).collect(Collectors.toList());

            List<Boolean> result = futures.stream().map(CompletableFuture::join).collect(Collectors.toList());

            if (result.contains(false)) {
                executor.shutdown();
                throw new IllegalStateException("There was an unhandled exception during chromosome interval import.");
            }
        }

        executor.shutdown();

        mDone = true;
    }

    private void updateConfigPartitionsAndLbUb(ImportConfig importConfig, final int index,
                                               final int rank) {
        String chromosomeName = importConfig.getImportConfiguration().getColumnPartitions(rank).getBegin()
                .getContigPosition().getContig();
        int chromosomeStart = (int) importConfig.getImportConfiguration().getColumnPartitions(rank).getBegin()
                .getContigPosition().getPosition();
        int chromosomeEnd = (int) importConfig.getImportConfiguration().getColumnPartitions(rank).getEnd()
                .getContigPosition().getPosition();
        GenomicsDBImportConfiguration.Partition partition = importConfig.getImportConfiguration()
            .getColumnPartitions(rank);
        if(partition.hasGenerateArrayNameFromPartitionBounds()) {
            String arrayName = String.format(CHROMOSOME_INTERVAL_FOLDER, chromosomeName, chromosomeStart, chromosomeEnd);
            partition = partition.toBuilder().setArrayName(arrayName).build();
        }
        GenomicsDBImportConfiguration.ImportConfiguration importConfiguration = importConfig
                .getImportConfiguration().toBuilder()
                .setLbCallsetRowIdx((long) index)
                .setUbCallsetRowIdx((long) (index + importConfig.getBatchSize() - 1))
                .setColumnPartitions(rank, partition).build();
        importConfig.setImportConfiguration(importConfiguration);
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
    public void write(final long lbRowIdx) throws GenomicsDBException {
        write(mLoaderJSONFile, mRank, lbRowIdx, Long.MAX_VALUE - 1);
    }

    /**
     * Write to TileDB/GenomicsDB using the configuration specified in the
     * loader file passed to constructor
     *
     * @param rank     Rank of this process (TileDB/GenomicsDB partition idx)
     * @param lbRowIdx Minimum row idx from which new data will be added
     */
    public void write(final int rank, final long lbRowIdx) throws GenomicsDBException {
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
    public void write(final int rank, final long lbRowIdx, final long ubRowIdx) throws GenomicsDBException {
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
    public void write(final String loaderJSONFile, final int rank, final long lbRowIdx, final long ubRowIdx)
            throws GenomicsDBException {
        mDone = false;
        if (loaderJSONFile == null) throw new GenomicsDBException("Loader JSON file not specified");
        if (mContainsBufferStreams) throw new GenomicsDBException("Cannot call write() functions if buffer streams are added");
        int status = jniGenomicsDBImporter(loaderJSONFile, rank, lbRowIdx, ubRowIdx);
        if (status != 0) throw new GenomicsDBException("GenomicsDBImporter write failed for loader JSON: "
                + loaderJSONFile + " rank: " + rank);
        mDone = true;
    }
}
