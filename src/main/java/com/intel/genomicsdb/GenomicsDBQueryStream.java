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

import java.io.InputStream;
import java.io.IOException;

/**
 * Provides a java.io.InputStream interface for the GenomicsDB combine gVCF operation.
 * It can be used as to construct a
 * <a href="https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/tribble/readers/PositionalBufferedStream.html">
 *   PositionalBufferedStream</a> object.The PositionalBufferedStream object can
 *   then be used by FeatureCodecs such as BCF2Codec to construct VariantContext objects
 */
public class GenomicsDBQueryStream extends InputStream
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
            throw new GenomicsDBException("Could not load genomicsdb native library", ule);
        }
    }
    private static boolean mDefaultReadAsBCF = true;
    private static boolean mDefaultUseMissingOnlyNotVectorEnd = true;
    private static boolean mDefaultKeepIDXFieldsInHeader = false;

    /*
     * Returns a "pointer" to a structure that stores the TileDB/GenomicsDB read state
     * This might look scary, but follows the same idea used in Java's compression library
     */
    private native long jniGenomicsDBInit(String loaderJSONFile, String queryJSONFile,
            String chr, int start, int end,
            int rank, long bufferCapacity, long segmentSize,
            boolean readAsBCF, boolean produceHeaderOnly,
            boolean useMissingValuesOnlyNotVectorEnd, boolean keepIDXFieldsInHeader);
    
    private native long jniGenomicsDBClose(long handle);

    private native long jniGenomicsDBGetNumBytesAvailable(long handle);

    private native byte jniGenomicsDBReadNextByte(long handle);

    private native int jniGenomicsDBRead(long handle, byte[] buffer, int off, int len);

    private native long jniGenomicsDBSkip(long handle, long n);

    private String mLoaderJSONFile;
    private String mQueryJSONFile;

    //"Pointer" to TileDB/GenomicsDB read state object
    private long mGenomicsDBReadStateHandle = 0;

    /**
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param queryJSONFile GenomicsDB query JSON configuration file
     */
    public GenomicsDBQueryStream(final String loaderJSONFile, final String queryJSONFile)
    {
        this(loaderJSONFile, queryJSONFile, "", 0, 0);
    }

    /**
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param queryJSONFile GenomicsDB query JSON configuration file
     * @param readAsBCF serialize-deserialize VCF records as BCF2 records
     */
    public GenomicsDBQueryStream(final String loaderJSONFile, final String queryJSONFile,
            final boolean readAsBCF)
    {
        this(loaderJSONFile, queryJSONFile, readAsBCF, false);
    }

    /**
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param queryJSONFile GenomicsDB query JSON configuration file
     * @param readAsBCF serialize-deserialize VCF records as BCF2 records
     * @param produceHeaderOnly produce only the header
     */
    public GenomicsDBQueryStream(final String loaderJSONFile, final String queryJSONFile,
            final boolean readAsBCF, final boolean produceHeaderOnly)
    {
        this(loaderJSONFile, queryJSONFile, "", 0, 0,
                0, 10485760,10485760,
                readAsBCF, produceHeaderOnly);
    }

    /**
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param queryJSONFile GenomicsDB query JSON configuration file
     * @param chr contig name
     * @param start start position (1-based)
     * @param end end position, inclusive (1-based)
     */
    public GenomicsDBQueryStream(final String loaderJSONFile, final String queryJSONFile,
            final String chr, final int start, final int end)
    {
        this(loaderJSONFile, queryJSONFile, chr, start, end, 0);
    }

    /**
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param queryJSONFile GenomicsDB query JSON configuration file
     * @param chr contig name
     * @param start start position (1-based)
     * @param end end position, inclusive (1-based)
     * @param readAsBCF serialize-deserialize VCF records as BCF2 records
     */
    public GenomicsDBQueryStream(final String loaderJSONFile, final String queryJSONFile,
            final String chr, final int start, final int end,
            final boolean readAsBCF)
    {
        this(loaderJSONFile, queryJSONFile, chr, start, end, 0, 10485760,
          10485760, readAsBCF, false);
    }
 
    /**
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param queryJSONFile GenomicsDB query JSON configuration file
     * @param chr contig name
     * @param start start position (1-based)
     * @param end end position, inclusive (1-based)
     * @param rank rank of this object if launched from within an MPI context (not used)
     */
    public GenomicsDBQueryStream(final String loaderJSONFile, final String queryJSONFile,
            final String chr, final int start, final int end,
            final int rank)
    {
        this(loaderJSONFile, queryJSONFile, chr, start, end, rank, 10485760,
          10485760);
    }

    /**
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param queryJSONFile GenomicsDB query JSON configuration file
     * @param chr contig name
     * @param start start position (1-based)
     * @param end end position, inclusive (1-based)
     * @param rank rank of this object if launched from within an MPI context (not used)
     * @param bufferCapacity size of buffer in bytes to be used by the native layer
     *                       to store combined BCF2 records
     * @param segmentSize buffer to be used for querying TileDB
     */
    public GenomicsDBQueryStream(final String loaderJSONFile, final String queryJSONFile,
            final String chr, final int start, final int end,
            final int rank, final long bufferCapacity, final long segmentSize)
    {
        this(loaderJSONFile, queryJSONFile, chr, start, end, rank, bufferCapacity,
          segmentSize, mDefaultReadAsBCF, false);
    }

    /**
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param queryJSONFile GenomicsDB query JSON configuration file
     * @param chr contig name
     * @param start start position (1-based)
     * @param end end position, inclusive (1-based)
     * @param rank rank of this object if launched from within an MPI context (not used)
     * @param bufferCapacity size of buffer in bytes to be used by the native layer
     *                       to store combined BCF2 records
     * @param segmentSize buffer to be used for querying TileDB
     * @param readAsBCF serialize-deserialize VCF records as BCF2 records
     * @param produceHeaderOnly produce VCF/BCF header only - no records (minor optimization)
     */
    public GenomicsDBQueryStream(final String loaderJSONFile, final String queryJSONFile,
            final String chr, final int start, final int end,
            final int rank, final long bufferCapacity, final long segmentSize,
            final boolean readAsBCF, final boolean produceHeaderOnly)
    {
        this(loaderJSONFile, queryJSONFile, chr, start, end, rank, bufferCapacity,
          segmentSize, readAsBCF, produceHeaderOnly,
          mDefaultUseMissingOnlyNotVectorEnd, mDefaultKeepIDXFieldsInHeader);
    }

    /**
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param queryJSONFile GenomicsDB query JSON configuration file
     * @param chr contig name
     * @param start start position (1-based)
     * @param end end position, inclusive (1-based)
     * @param rank rank of this object if launched from within an MPI context (not used)
     * @param bufferCapacity size of buffer in bytes to be used by the native layer
     *                       to store combined BCF2 records
     * @param segmentSize buffer to be used for querying TileDB
     * @param readAsBCF serialize-deserialize VCF records as BCF2 records
     * @param produceHeaderOnly produce VCF/BCF header only - no records (minor optimization)
     * @param useMissingValuesOnlyNotVectorEnd don't add BCF2.2 vector end values
     * @param keepIDXFieldsInHeader keep BCF IDX fields in header
     */
    public GenomicsDBQueryStream(final String loaderJSONFile, final String queryJSONFile,
            final String chr, final int start, final int end,
            final int rank, final long bufferCapacity, final long segmentSize,
            final boolean readAsBCF, final boolean produceHeaderOnly,
            final boolean useMissingValuesOnlyNotVectorEnd, final boolean keepIDXFieldsInHeader)
    {
        mLoaderJSONFile = loaderJSONFile;
        mQueryJSONFile = queryJSONFile;
        mGenomicsDBReadStateHandle = jniGenomicsDBInit(loaderJSONFile, queryJSONFile, chr,
          start, end, rank, bufferCapacity, segmentSize,
          readAsBCF, produceHeaderOnly,
          useMissingValuesOnlyNotVectorEnd, keepIDXFieldsInHeader);
    }

    @Override
    public int available() throws IOException
    {
        return (int)jniGenomicsDBGetNumBytesAvailable(mGenomicsDBReadStateHandle);
    }

    @Override
    public void close() throws IOException
    {
        mGenomicsDBReadStateHandle = jniGenomicsDBClose(mGenomicsDBReadStateHandle);
    }

    @Override
    public boolean markSupported()
    {
        return false;
    }

    @Override
    public int read() throws IOException
    {
        return jniGenomicsDBReadNextByte(mGenomicsDBReadStateHandle);
    }

    @Override
    public int read(byte[] buffer, int off, int len) throws IOException
    {
        if(len <= 0)
            return 0;
        long numBytesRead = jniGenomicsDBRead(mGenomicsDBReadStateHandle, buffer, off, len);
        return (numBytesRead == 0) ? -1 : (int)numBytesRead;
    }

    @Override
    public long skip(long n) throws IOException
    {
        return jniGenomicsDBSkip(mGenomicsDBReadStateHandle, n);
    }
}
