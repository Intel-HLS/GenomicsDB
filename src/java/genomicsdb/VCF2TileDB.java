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

    /**
     * JNI function that writes to TileDB
     */
    private native int jniVCF2TileDB(String loaderJSONFile, int rank, long lbRowIdx, long ubRowIdx);

    private String mLoaderJSONFile = null;
    private int mRank = 0;

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
        mLoaderJSONFile = loaderJSONFile;
    }

    /** 
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param rank Rank of this process (TileDB/GenomicsDB partition idx)
     */
    public VCF2TileDB(String loaderJSONFile, int rank)
    {
        mLoaderJSONFile = loaderJSONFile;
        mRank = rank;
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
        if(loaderJSONFile == null)
            throw new GenomicsDBException("Loader JSON file not specified");
        int status = jniVCF2TileDB(loaderJSONFile, rank, lbRowIdx, ubRowIdx);
        if(status != 0)
            throw new GenomicsDBException("VCF2TileDB write failed for loader JSON: "+loaderJSONFile+" rank: "+rank);
    }
}
