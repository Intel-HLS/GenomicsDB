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

package genomicsdb;

import java.lang.Long;
import java.lang.Integer;
import java.io.OutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.List;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.samtools.util.RuntimeIOException;
import org.json.simple.JSONObject;
import org.json.simple.JSONArray;
import java.io.StringWriter;

/**
 * Java wrapper for vcf2tiledb - imports VCFs into TileDB/GenomicsDB.
 * All vid information is assumed to be set correctly by the user (JSON files)
 */

public class VCF2TileDB
{
    static
    {
        try
        {
            boolean loaded = GenomicsDBUtils.loadLibrary();
            if(!loaded)
                throw new GenomicsDBException("Could not load genomicsdb native library");
        }
        catch(UnsatisfiedLinkError ule)
        {
            throw new GenomicsDBException("Could not load genomicsdb native library");
        }
    }
    private static long mDefaultBufferCapacity = 20480; //20KB

    /**
     * Buffer stream implementation - it's silent in the sense that when the buffer is full, it doesn't raise an exception
     * but just marks a a flag as full. It's up to the caller to check the flag and retry later
     * Why? Most likely, it's faster to check a flag rather than throw and catch an exception
     */
    private class SilentByteBufferStream extends OutputStream
    {
        private byte mBuffer[] = null;
        private long mNumValidBytes = 0;
        private long mMarker = 0;
        private boolean mOverflow = false;
        
        public SilentByteBufferStream()
        {
            mBuffer = new byte[(int)mDefaultBufferCapacity];
        }

        public SilentByteBufferStream(final long capacity)
        {
            mBuffer = new byte[(int)capacity];
        }

        @Override
        public void close() throws IOException  //does nothing
        {
        }

        @Override
        public void flush() throws IOException  //does nothing
        {
        }

        @Override
        public void write(byte[] b, int off, int len) throws IOException
        {
            if(mOverflow)
                return;
            if(len+mNumValidBytes > mBuffer.length)
                mOverflow = true;
            else
            {
                System.arraycopy(b, off, mBuffer, (int)mNumValidBytes, len);
                mNumValidBytes += len;
            }
        }

        @Override
        public void write(byte[] b) throws IOException
        {
            write(b, 0, b.length);
        }

        @Override
        public void write(int b) throws IOException
        {
            if(mOverflow)
                return;
            if(mNumValidBytes+1>mBuffer.length)
                mOverflow = true;
            else
            {
                mBuffer[(int)mNumValidBytes] = (byte)b;
                ++mNumValidBytes;
            }
        }

        public int size()
        {
            return mBuffer.length;
        }

        public void resize(final long newSize)
        {
            byte tmp[] = new byte[(int)newSize];
            System.arraycopy(mBuffer, 0, tmp, 0, mBuffer.length);
            mBuffer = tmp; //hopefully Java GC does its job
        }

        public boolean overflow()
        {
            return mOverflow;
        }

        public void setOverflow(final boolean value)
        {
            mOverflow = value;
        }

        public long getNumValidBytes()
        {
            return mNumValidBytes;
        }

        public void setNumValidBytes(final long value)
        {
            mNumValidBytes = value;
        }

        public void setMarker(final long value)
        {
            mMarker = value;
        }

        public long getMarker()
        {
            return mMarker;
        }

        public byte[] getBuffer()
        {
            return mBuffer;
        }
    }

    /**
     * Utility class wrapping a stream and a VariantContextWriter for a given stream
     * Each GenomicsDB import stream consists of a buffer stream and a writer object
     */
    private class GenomicsDBImporterStreamWrapper
    {
        public VariantContextWriter mVCWriter = null;
        public SilentByteBufferStream mStream = null;
        private Iterator<VariantContext> mIterator = null;
        private VariantContext mCurrentVC = null;
        
        /**
         * Constructor
         * @param vcfHeader VCF header for the stream
         * @param bufferCapacity Capacity of the stream buffer in bytes
         * @param streamType BCF_STREAM or VCF_STREAM
         */
        public GenomicsDBImporterStreamWrapper(final VCFHeader vcfHeader, final long bufferCapacity,
                final VariantContextWriterBuilder.OutputType streamType, Iterator<VariantContext> vcIterator) throws GenomicsDBException
        {
            mIterator = vcIterator;
            if(vcIterator != null && vcIterator.hasNext())
                mCurrentVC = vcIterator.next();
            boolean headerWritten = false;
            long currentCapacity = bufferCapacity;
            //Must ensure that the header gets written into the buffer stream
            //Why this big outer loop? VCFWriter/BCFWriter seems to store some state which makes
            //calling writeHeader() multiple times impossible
            //Hence, create new objects in every iteration of the loop
            //Since this function is called only once per stream, not really
            //a performance concern
            while(!headerWritten)
            {
                mStream = new SilentByteBufferStream(currentCapacity);
                switch(streamType)
                {
                    case BCF_STREAM:
                        mVCWriter = new VariantContextWriterBuilder().setOutputBCFStream(mStream).unsetOption(Options.INDEX_ON_THE_FLY).build();
                        break;
                    case VCF_STREAM:
                        mVCWriter = new VariantContextWriterBuilder().setOutputVCFStream(mStream).unsetOption(Options.INDEX_ON_THE_FLY).build();
                        break;
                    default:
                        throw new GenomicsDBException("Unknown stream type "+streamType.toString());
                }
                //Why clone the header?
                //The writer modifies the VCFHeader object passed to writeHeader() - however, we might need
                //to call writeHeader multiple times if the underlying buffer in mStream is too small. Hence,
                //always pass a clone of the original, unmodified header in each call of writeHeader
                mVCWriter.writeHeader(new VCFHeader(vcfHeader));
                if(mStream.overflow())
                    currentCapacity = 2*currentCapacity+1;
                else
                    headerWritten = true;
            }
        }

        public boolean hasIterator()
        {
            return (mIterator != null);
        }

        public VariantContext next()
        {
            if(mIterator != null && mIterator.hasNext())
                mCurrentVC = mIterator.next();
            else
                mCurrentVC = null;
            return mCurrentVC;
        }

        public VariantContext getCurrentVC()
        {
            return mCurrentVC;
        }
    }

    /**
     * Utility class that stores row index and globally unique name for a given sample
     */
    public static class SampleInfo
    {
        public String mName = null;
        public long mRowIdx = -1;

        public SampleInfo(final String name, final long rowIdx)
        {
            mName = name;
            mRowIdx = rowIdx;
        }
    }

    /**
     * JNI functions
     */
    /**
     * Creates VCF2TileDB object when importing VCF files (no streams)
     * @param loaderJSONFile Path to loader JSON file
     * @param rank Rank of object - corresponds to the partition index in the loader for which this object will import data
     * @param lbRowIdx Smallest row idx which should be imported by this object
     * @param ubRowIdx Largest row idx which should be imported by this object
     * @return status - 0 if everything was ok, -1 otherwise
     */
    private native int jniVCF2TileDB(String loaderJSONFile, int rank, long lbRowIdx, long ubRowIdx);
    /**
     * Creates VCF2TileDB object when importing VCF files (no streams)
     * @param loaderJSONFile Path to loader JSON file
     * @param rank Rank of object - corresponds to the partition index in the loader for which this object will import data
     * @param lbRowIdx Smallest row idx which should be imported by this object
     * @param ubRowIdx Largest row idx which should be imported by this object
     * @return "pointer"/address to GenomicsDBImporter object in memory, if 0, then something went wrong
     */
    private native long jniInitializeGenomicsDBImporterObject(String loaderJSONFile, int rank, long lbRowIdx, long ubRowIdx);
    /**
     * Notify importer object that a new stream is to be added
     * @param genomicsDBImporterHandle "pointer" returned by jniInitializeGenomicsDBImporterObject
     * @param isBCF
     * @param bufferCapacity in bytes
     * @param buffer initialization buffer containing the VCF/BCF header
     * @param numValidBytesInBuffer num valid bytes in the buffer (length of the header)
     */
    private native void jniAddBufferStream(long genomicsDBImporterHandle, String streamName, boolean isBCF,
            long bufferCapacity, byte[] buffer, long numValidBytesInBuffer);
    /**
     * Setup loader after all the buffer streams are added
     * @param genomicsDBImporterHandle "pointer" returned by jniInitializeGenomicsDBImporterObject
     * @param callsetMappingJSON JSON formatted string containing globally consistent callset name to row index mapping
     * @return maximum number of buffer stream identifiers that can be returned in mExhaustedBufferStreamIdentifiers later
     * (this depends on the number of partitions and the number of buffer streams)
     */
    private native long jniSetupGenomicsDBLoader(long genomicsDBImporterHandle, final String callsetMappingJSON);
    /**
     * @param genomicsDBImporterHandle "pointer" returned by jniInitializeGenomicsDBImporterObject
     * @param streamIdx stream index
     * @param partitionIdx partition index (unused now)
     * @param buffer buffer containing data
     * @param numValidBytesInBuffer num valid bytes in the buffer
     */
    private native void jniWriteDataToBufferStream(long handle, int streamIdx, int partitionIdx,
            byte[] buffer, long numValidBytesInBuffer);
    /**
     * Import the next batch of data into TileDB/GenomicsDB
     * @param genomicsDBImporterHandle "pointer" returned by jniInitializeGenomicsDBImporterObject
     * @param exhaustedBufferIdentifiers contains the list of exhausted buffer stream identifiers - the number of 
     * exhausted streams is stored in the last element of the array
     * @return true if the whole import process is completed, false otherwise
     */
    private native boolean jniImportBatch(long genomicsDBImporterHandle, long[] exhaustedBufferIdentifiers);

    private String mLoaderJSONFile = null;
    private int mRank = 0;
    private long mLbRowIdx = 0;
    private long mUbRowIdx = Long.MAX_VALUE-1;

    //For buffered streams
    private boolean mContainsBufferStreams = false;
    private long mGenomicsDBImporterObjectHandle = 0;
    private ArrayList<GenomicsDBImporterStreamWrapper> mBufferStreamWrapperVector = null;
    private boolean mIsLoaderSetupDone = false;
    //To find out which buffer streams are exhausted
    private long mMaxBufferStreamIdentifiers = 0;
    private long[] mExhaustedBufferStreamIdentifiers = null;
    private long mNumExhaustedBufferStreams = 0;
    //Done flag - useful only for buffered streams
    private boolean mDone = false;
    //JSON object that specifies callset/sample name to row_idx mapping in the buffer
    private JSONObject mCallsetMappingJSON = null;

    /**
     * Constructor
     */
    public VCF2TileDB()
    {
    }

    /** 
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     */
    public VCF2TileDB(String loaderJSONFile)
    {
        initialize(loaderJSONFile, 0, 0, Long.MAX_VALUE-1);
    }

    /** 
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank Rank of this process (TileDB/GenomicsDB partition idx)
     */
    public VCF2TileDB(String loaderJSONFile, int rank)
    {
        initialize(loaderJSONFile, rank, 0, Long.MAX_VALUE-1);
    }
    /** 
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank Rank of this process (TileDB/GenomicsDB partition idx)
     * @param lbRowIdx Smallest row idx which should be imported by this object
     * @param ubRowIdx Largest row idx which should be imported by this object
     */
    public VCF2TileDB(String loaderJSONFile, int rank, long lbRowIdx, long ubRowIdx)
    {
        initialize(loaderJSONFile, rank, lbRowIdx, ubRowIdx);
    }

    /**
     * Initialize variables
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank Rank of this process (TileDB/GenomicsDB partition idx)
     * @param lbRowIdx Smallest row idx which should be imported by this object
     * @param ubRowIdx Largest row idx which should be imported by this object
     */
    private void initialize(String loaderJSONFile, int rank, long lbRowIdx, long ubRowIdx)
    {
        mLoaderJSONFile = loaderJSONFile;
        mRank = rank;
        mLbRowIdx = lbRowIdx;
        mUbRowIdx = ubRowIdx;
    }

    /**
     * Static function that reads sample names from the vcfHeader and adds entries to the map.The function assumes that
     * the samples will be assigned row indexes beginning at rowIdx and that the sample names specified in the header
     * are globally unique (across all streams/files)
     * @param sampleIndexToInfo  <sampleIndex in vcfHeader: SampleInfo> map
     * @param vcfHeader VCF header
     * @param rowIdx Starting row index from which to assign
     * @return rowIdx+#samples in the header
     */
    public static long initializeSampleInfoMapFromHeader(Map<Integer, SampleInfo> sampleIndexToInfo, final VCFHeader vcfHeader, final long rowIdx)
    {
        final List<String> headerSampleNames = vcfHeader.getGenotypeSamples();
        final int numSamplesInHeader = headerSampleNames.size();
        for(int i=0;i<numSamplesInHeader;++i)
            sampleIndexToInfo.put(i, new SampleInfo(headerSampleNames.get(i), rowIdx+i));
        return rowIdx + numSamplesInHeader;
    }

    /**
     * Add a buffer stream as the data source - caller must:
     * 1. Call setupGenomicsDBImporter() once all streams are added
     * 2. Provide VC objects using the add() function
     * 3. Call importBatch()
     * 4. Get list of exhausted buffer streams using getNumExhaustedBufferStreams(), getExhaustedBufferStreamIndex()
     * 5. If !isDone() goto 2
     * @param streamName Name of the stream being added - must be unique with respect to this VCF2TileDB object
     * @param vcfHeader VCF header for the stream
     * @param bufferCapacity Capacity of the stream buffer in bytes
     * @param streamType BCF_STREAM or VCF_STREAM
     * @param sampleIndexToInfo map from sample index in the vcfHeader to SampleInfo object which contains row index and globally unique name
     * can be set to null, which implies that the mapping is stored in a callsets JSON file
     */
    public int addBufferStream(final String streamName, final VCFHeader vcfHeader, final long bufferCapacity,
            final VariantContextWriterBuilder.OutputType streamType,
            final Map<Integer, SampleInfo> sampleIndexToInfo) throws GenomicsDBException
    {
        return addBufferStream(streamName, vcfHeader, bufferCapacity, streamType, null, sampleIndexToInfo);
    }

    /**
     * Add a sorted VC iterator as the data source - caller must:
     * 1. Call setupGenomicsDBImporter() once all iterators are added
     * 2. Call importBatch()
     * 3. Done!
     * @param streamName Name of the stream being added - must be unique with respect to this VCF2TileDB object
     * @param vcfHeader VCF header for the stream
     * @param vcIterator Iterator over VariantContext objects
     * @param bufferCapacity Capacity of the stream buffer in bytes
     * @param streamType BCF_STREAM or VCF_STREAM
     * @param sampleIndexToInfo map from sample index in the vcfHeader to SampleInfo object which contains row index and globally unique name
     * can be set to null, which implies that the mapping is stored in a callsets JSON file
     */
    public int addSortedVariantContextIterator(final String streamName, final VCFHeader vcfHeader, Iterator<VariantContext> vcIterator,
            final long bufferCapacity, final VariantContextWriterBuilder.OutputType streamType,
            final Map<Integer, SampleInfo> sampleIndexToInfo) throws GenomicsDBException
    {
        return addBufferStream(streamName, vcfHeader, bufferCapacity, streamType, vcIterator, sampleIndexToInfo);
    }

    /**
     * Sets sorted VC iterator as the data source and calls setupGenomicsDBImporter(). No more streams/iterators can
     * be added after this function is called. The caller must:
     * 1. Call importBatch()
     * 2. Done!
     * @param streamName Name of the stream being added - must be unique with respect to this VCF2TileDB object
     * @param vcfHeader VCF header for the stream
     * @param vcIterator Iterator over VariantContext objects
     * @param bufferCapacity Capacity of the stream buffer in bytes
     * @param streamType BCF_STREAM or VCF_STREAM
     * @param sampleIndexToInfo map from sample index in the vcfHeader to SampleInfo object which contains row index and globally unique name
     * can be set to null, which implies that the mapping is stored in a callsets JSON file
     */
    public int setSortedVariantContextIterator(final String streamName, final VCFHeader vcfHeader, Iterator<VariantContext> vcIterator,
            final long bufferCapacity, final VariantContextWriterBuilder.OutputType streamType,
            final Map<Integer, SampleInfo> sampleIndexToInfo) throws GenomicsDBException, IOException
    {
        int streamIdx = addSortedVariantContextIterator(streamName, vcfHeader, vcIterator, bufferCapacity, streamType, sampleIndexToInfo);
        setupGenomicsDBImporter();
        return streamIdx;
    }

    /**
     * Add a buffer stream or VC iterator - internal function
     * @param streamName Name of the stream being added - must be unique with respect to this VCF2TileDB object
     * @param vcfHeader VCF header for the stream
     * @param bufferCapacity Capacity of the stream buffer in bytes
     * @param streamType BCF_STREAM or VCF_STREAM
     * @param vcIterator Iterator over VariantContext objects - can be null
     * @param sampleIndexToInfo map from sample index in the vcfHeader to SampleInfo object which contains row index and globally unique name
     * can be set to null, which implies that the mapping is stored in a callsets JSON file
     */
    private int addBufferStream(final String streamName, final VCFHeader vcfHeader, final long bufferCapacity,
            final VariantContextWriterBuilder.OutputType streamType, Iterator<VariantContext> vcIterator,
            final Map<Integer, SampleInfo> sampleIndexToInfo) throws GenomicsDBException
    {
        if(mIsLoaderSetupDone)
            throw new GenomicsDBException("Cannot add buffer streams after setupGenomicsDBImporter() is called");
        //First time a buffer is added
        if(!mContainsBufferStreams)
        {
            mGenomicsDBImporterObjectHandle = jniInitializeGenomicsDBImporterObject(mLoaderJSONFile, mRank, mLbRowIdx, mUbRowIdx);
            if(mGenomicsDBImporterObjectHandle == 0)
                throw new GenomicsDBException("Could not initialize GenomicsDBImporter object");
            mBufferStreamWrapperVector = new ArrayList<GenomicsDBImporterStreamWrapper>();
            mCallsetMappingJSON = new JSONObject();
            mContainsBufferStreams = true;
        }
        mBufferStreamWrapperVector.add(new GenomicsDBImporterStreamWrapper(vcfHeader, bufferCapacity, streamType, vcIterator));
        int currIdx = mBufferStreamWrapperVector.size()-1;
        SilentByteBufferStream currStream = mBufferStreamWrapperVector.get(currIdx).mStream;
        jniAddBufferStream(mGenomicsDBImporterObjectHandle, streamName, streamType == VariantContextWriterBuilder.OutputType.BCF_STREAM,
                bufferCapacity, currStream.getBuffer(), currStream.getNumValidBytes());
        if(sampleIndexToInfo != null)
        {
            for(Map.Entry<Integer, SampleInfo> currEntry : sampleIndexToInfo.entrySet())
            {
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
     * Setup the importer after all the buffer streams are added, but before any data is inserted into any stream
     * No more buffer streams can be added once setupGenomicsDBImporter() is called
     */
    public void setupGenomicsDBImporter() throws IOException
    {
        if(mIsLoaderSetupDone)
            return;
        //Callset mapping JSON - convert to string
        JSONObject topCallsetJSON = new JSONObject();
        topCallsetJSON.put("callsets", mCallsetMappingJSON);
        StringWriter stringWriter = new StringWriter();
        topCallsetJSON.writeJSONString(stringWriter);
        //Call native setupGenomicsDBImporter()
        mMaxBufferStreamIdentifiers = jniSetupGenomicsDBLoader(mGenomicsDBImporterObjectHandle, stringWriter.toString());
        //Why 2* - each identifier is a pair<buffer_stream_idx, partition_idx>
        //Why +1 - the last element will contain the number of exhausted stream identifiers when importBatch() is called
        mExhaustedBufferStreamIdentifiers = new long[2*((int)mMaxBufferStreamIdentifiers)+1];
        //Set all streams to empty
        //Add all streams to mExhaustedBufferStreamIdentifiers - this way when importBatch is called the first
        //time all streams' data are written
        for(int i=0,idx=0;i<mBufferStreamWrapperVector.size();++i,idx+=2)
        {
            SilentByteBufferStream currStream = mBufferStreamWrapperVector.get(i).mStream;
            currStream.setNumValidBytes(0);
            mExhaustedBufferStreamIdentifiers[idx] = i;
            mExhaustedBufferStreamIdentifiers[idx+1] = 0;
        }
        //Set number of exhausted buffer streams - all streams are exhausted the first time
        mNumExhaustedBufferStreams = mBufferStreamWrapperVector.size();
        mExhaustedBufferStreamIdentifiers[(int)(2*mMaxBufferStreamIdentifiers)] = mNumExhaustedBufferStreams;
        mIsLoaderSetupDone = true;
    }

    /**
     * Write VariantContext object to stream - may fail if the buffer is full
     * It's the caller's responsibility keep track of the VC object that's not written
     * @param vc VariantContext object
     * @param streamIdx index of the stream returned by the addBufferStream() call
     * @return true if the vc object was written successfully, false otherwise
     */
    public boolean add(VariantContext vc, final int streamIdx) throws GenomicsDBException, RuntimeIOException
    {
        if(streamIdx < 0 || streamIdx >= mBufferStreamWrapperVector.size())
            throw new GenomicsDBException("Invalid stream idx "+Integer.toString(streamIdx)+" must be between [0-"
                    +Long.toString(mBufferStreamWrapperVector.size()-1)+"]");
        if(!mIsLoaderSetupDone)
            throw new GenomicsDBException("Cannot add VariantContext objects to streams before calling setupGenomicsDBImporter()");
        GenomicsDBImporterStreamWrapper currWrapper = mBufferStreamWrapperVector.get(streamIdx);
        currWrapper.mVCWriter.add(vc);
        SilentByteBufferStream currStream = currWrapper.mStream;
        if(currStream.overflow() && currStream.getMarker() > 0) //at least one record already existed in the buffer
        {
            //Set num valid bytes to marker - marker points to location after the last valid serialized vc in the buffer
            currStream.setNumValidBytes(currStream.getMarker());
            return false;
        }
        else
        {
            while(currStream.overflow()) //the first record to be added to the buffer is too large, resize buffer
            {
                currStream.resize(2*currStream.size()+1);
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
     */
    public boolean importBatch() throws IOException
    {
        if(mDone)
            return true;
        if(!mIsLoaderSetupDone)
            setupGenomicsDBImporter();
        boolean allExhaustedStreamsHaveIterators = true;
        while(!mDone && allExhaustedStreamsHaveIterators)
        {
            //Write data from buffer streams exhausted in the previous round into GenomicsDB
            for(int i=0,idx=0;i<mNumExhaustedBufferStreams;++i,idx+=2)
            {
                int bufferStreamIdx = (int)mExhaustedBufferStreamIdentifiers[idx];
                GenomicsDBImporterStreamWrapper currWrapper = mBufferStreamWrapperVector.get(bufferStreamIdx);
                //If iterator is provided, get data from iterator
                if(currWrapper.hasIterator())
                {
                    while(currWrapper.getCurrentVC() != null)
                    {
                        boolean added = add(currWrapper.getCurrentVC(), bufferStreamIdx);
                        if(added)
                            currWrapper.next();
                        else
                            break; //buffer full
                    }
                }
                SilentByteBufferStream currStream = currWrapper.mStream;
                jniWriteDataToBufferStream(mGenomicsDBImporterObjectHandle, bufferStreamIdx, 0, currStream.getBuffer(), currStream.getNumValidBytes());
            }
            mDone = jniImportBatch(mGenomicsDBImporterObjectHandle, mExhaustedBufferStreamIdentifiers);
            mNumExhaustedBufferStreams = mExhaustedBufferStreamIdentifiers[mExhaustedBufferStreamIdentifiers.length-1];
            //Reset markers, numValidBytesInBuffer and overflow flag for the exhausted streams
            for(long i=0, idx=0;i<mNumExhaustedBufferStreams;++i, idx+=2)
            {
                int bufferStreamIdx = (int)mExhaustedBufferStreamIdentifiers[(int)idx];
                GenomicsDBImporterStreamWrapper currWrapper = mBufferStreamWrapperVector.get(bufferStreamIdx);
                if(!currWrapper.hasIterator())
                    allExhaustedStreamsHaveIterators = false;
                SilentByteBufferStream currStream = currWrapper.mStream;
                currStream.setOverflow(false);
                currStream.setMarker(0);
                currStream.setNumValidBytes(0);
            }
            if(mDone)
            {
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
    public long getNumExhaustedBufferStreams()
    {
        return mNumExhaustedBufferStreams;
    }

    /**
     * Get buffer stream index of i-th exhausted stream
     * There are mNumExhaustedBufferStreams and the caller must provide data for streams
     * with indexes getExhaustedBufferStreamIndex(0), getExhaustedBufferStreamIndex(1),..., getExhaustedBufferStreamIndex(mNumExhaustedBufferStreams-1)
     * @param i i-th exhausted buffer stream
     */
    public int getExhaustedBufferStreamIndex(final long i)
    {
        assert i < mNumExhaustedBufferStreams && i >= 0;
        return (int)mExhaustedBufferStreamIdentifiers[2*((int)i)]; //why 2* - exhausted buffer stream identifier is a pair<stream_idx, partition_idx>
    }

    /**
     * Is the import process completed
     * @return true if complete, false otherwise
     */
    public boolean isDone()
    {
        return mDone;
    }

    /**
     * Write to TileDB/GenomicsDB using the configuration specified in the
     * loader file passed to constructor
     */
    public void write() throws GenomicsDBException
    {
        write(mLoaderJSONFile, mRank, 0, Long.MAX_VALUE-1);
    }

    /**
     * Write to TileDB/GenomicsDB using the configuration specified in the
     * loader file passed to constructor
     * @param lbRowIdx Minimum row idx from which new data will be added
     */
    public void write(long lbRowIdx) throws GenomicsDBException
    {
        write(mLoaderJSONFile, mRank, lbRowIdx, Long.MAX_VALUE-1);
    }

    /**
     * Write to TileDB/GenomicsDB using the configuration specified in the
     * loader file passed to constructor
     * @param rank Rank of this process (TileDB/GenomicsDB partition idx)
     * @param lbRowIdx Minimum row idx from which new data will be added
     */
    public void write(int rank, long lbRowIdx) throws GenomicsDBException
    {
        write(mLoaderJSONFile, rank, lbRowIdx, Long.MAX_VALUE-1);
    }

    /**
     * Write to TileDB/GenomicsDB using the configuration specified in the
     * loader file passed to constructor
     * @param rank Rank of this process (TileDB/GenomicsDB partition idx)
     * @param lbRowIdx Minimum row idx from which new data will be added
     * @param ubRowIdx Maximum row idx upto which new data will be added
     */
    public void write(int rank, long lbRowIdx, long ubRowIdx) throws GenomicsDBException
    {
        write(mLoaderJSONFile, rank, lbRowIdx, ubRowIdx);
    }

    /**
     * Write to TileDB/GenomicsDB
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank Rank of this process (TileDB/GenomicsDB partition idx)
     * @param lbRowIdx Minimum row idx from which new data will be added
     * @param ubRowIdx Maximum row idx upto which new data will be added
     */
    public void write(String loaderJSONFile, int rank, long lbRowIdx, long ubRowIdx) throws GenomicsDBException
    {
        mDone = false;
        if(loaderJSONFile == null)
            throw new GenomicsDBException("Loader JSON file not specified");
        if(mContainsBufferStreams)
            throw new GenomicsDBException("Cannot call write() functions if buffer streams are added");
        int status = jniVCF2TileDB(loaderJSONFile, rank, lbRowIdx, ubRowIdx);
        if(status != 0)
            throw new GenomicsDBException("VCF2TileDB write failed for loader JSON: "+loaderJSONFile+" rank: "+rank);
        mDone = true;
    }
}
