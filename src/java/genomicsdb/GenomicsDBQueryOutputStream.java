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

import java.io.InputStream;
import java.io.IOException;

/**
 * Provides a java.io.InputStream interface for a JNI buffer.This class
 * can be used as to construct a <a href="https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/tribble/readers/PositionalBufferedStream.html">PositionalBufferedStream</a> object.The PositionalBufferedStream object can then be used by FeatureCodecs such as BCF2Codec to construct
 * VariantContext objects
 */
public class GenomicsDBQueryOutputStream extends InputStream
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
    /*
     * Returns a "pointer" to a structure that stores the TileDB/GenomicsDB read state
     * This might look scary, but follows the same idea used in Java's compression library
     */
    private native long jniGenomicsDBInit(String loaderJSONFile, String queryJSONFile,
            String chr, int start, int end,
            int rank, long bufferCapacity, long segmentSize);
    
    private native long jniGenomicsDBClose(long handle);

    private native long jniGenomicsDBGetNumBytesAvailable(long handle);

    private native byte jniGenomicsDBReadNextByte(long handle);

    private native int jniGenomicsDBRead(long handle, byte[] buffer, int off, int len);

    private native long jniGenomicsDBSkip(long handle, long n);

    private String mLoaderJSONFile;
    private String mQueryJSONFile;

    //"Pointer" to TileDB/GenomicsDB read state object
    private long mGenomicsDBReadStateHandle = 0;

    public GenomicsDBQueryOutputStream(final String loaderJSONFile, final String queryJSONFile)
    {
        this(loaderJSONFile, queryJSONFile, "", 0, 0);
    }

    public GenomicsDBQueryOutputStream(final String loaderJSONFile, final String queryJSONFile,
            final String chr, final int start, final int end)
    {
        this(loaderJSONFile, queryJSONFile, chr, start, end, 0);
    }

    public GenomicsDBQueryOutputStream(final String loaderJSONFile, final String queryJSONFile,
            final String chr, final int start, final int end,
            final int rank)
    {
        this(loaderJSONFile, queryJSONFile, chr, start, end, rank, 10485760, 10485760);
    }

    public GenomicsDBQueryOutputStream(final String loaderJSONFile, final String queryJSONFile,
            final String chr, final int start, final int end,
            final int rank, final long bufferCapacity, final long segmentSize)
    {
        mLoaderJSONFile = loaderJSONFile;
        mQueryJSONFile = queryJSONFile;
        mGenomicsDBReadStateHandle = jniGenomicsDBInit(loaderJSONFile, queryJSONFile, chr, start, end, rank, bufferCapacity, segmentSize);
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
