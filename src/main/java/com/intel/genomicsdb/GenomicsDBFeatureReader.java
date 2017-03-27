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

package com.intel.genomicsdb;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.bcf2.BCF2Codec;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * A reader for GenomicsDB that implements {@link htsjdk.tribble.FeatureReader}
 * Currently, the reader only return {@link htsjdk.variant.variantcontext.VariantContext}
 */

public class GenomicsDBFeatureReader<T extends Feature, SOURCE> implements FeatureReader<T>
{
    private String mLoaderJSONFile = null;
    private String mQueryJSONFile = null;
    private FeatureCodec<T, SOURCE> mCodec = null;
    protected FeatureCodecHeader mFCHeader;
    private VCFHeader mVCFHeader = null;
    private ArrayList<String> mSequenceNames;

    /**
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param tiledbWorkspace TileDB workspace path
     * @param arrayName TileDB array name
     * @param referenceGenome Path to reference genome (fasta file)
     * @param codec FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *              and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @throws IOException when data cannot be read from the stream 
     */
    public GenomicsDBFeatureReader(final String loaderJSONFile,
            final String tiledbWorkspace, final String arrayName,
            final String referenceGenome,
            final FeatureCodec<T, SOURCE> codec) throws IOException
    {
        this(loaderJSONFile, tiledbWorkspace, arrayName,
                referenceGenome, null,
                codec);
    }

    /**
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param tiledbWorkspace TileDB workspace path
     * @param arrayName TileDB array name
     * @param referenceGenome Path to reference genome (fasta file)
     * @param templateVCFHeaderFilename Template VCF header to be used for
     *                                  the combined gVCF records
     * @param codec FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *              and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @throws IOException when data cannot be read from the stream 
     */
    public GenomicsDBFeatureReader(final String loaderJSONFile,
            final String tiledbWorkspace, final String arrayName,
            final String referenceGenome, final String templateVCFHeaderFilename,
            final FeatureCodec<T, SOURCE> codec) throws IOException
    {
        //Produce temporary JSON query config file
        String indentString = "    ";
        String queryJSON = "{\n";
        queryJSON += indentString + "\"scan_full\": true,\n";
        queryJSON += indentString + "\"workspace\": \""+tiledbWorkspace+"\",\n";
        queryJSON += indentString + "\"array\": \""+arrayName+"\",\n";
        queryJSON += indentString + "\"reference_genome\": \""+referenceGenome+"\",\n";
        if(templateVCFHeaderFilename != null)
            queryJSON += indentString + "\"vcf_header_filename\": \"" +
              templateVCFHeaderFilename+"\"\n";
        queryJSON += "}\n";
        File tmpQueryJSONFile = File.createTempFile("queryJSON", ".json");
        tmpQueryJSONFile.deleteOnExit();
        FileWriter fptr = new FileWriter(tmpQueryJSONFile);
        fptr.write(queryJSON);
        fptr.close();
        initialize(loaderJSONFile, tmpQueryJSONFile.getAbsolutePath(), codec);
    }

    /**
     * Constructor
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param queryJSONFile GenomicsDB query JSON configuration file
     * @param codec FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *              and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @throws IOException when data cannot be read from the stream 
     */
    public GenomicsDBFeatureReader(final String loaderJSONFile, final String queryJSONFile,
            final FeatureCodec<T, SOURCE> codec) throws IOException
    {
        initialize(loaderJSONFile, queryJSONFile, codec);
    }

    /**
     * Constructor
     * @param queryJSONFile GenomicsDB query JSON configuration file. Since the constructor
     * has no loader JSON as argument, you must specify the vid and callset mapping in the query JSON
     * @param codec FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *              and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @throws IOException when data cannot be read from the stream
     */
    public GenomicsDBFeatureReader(final String queryJSONFile,
            final FeatureCodec<T, SOURCE> codec) throws IOException

    {
        initialize("", queryJSONFile, codec);
    }

    /**
     * Initialization function that's used by all constructors
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param queryJSONFile GenomicsDB query JSON configuration file
     * @param codec FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *              and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @throws IOException when data cannot be read from the stream 
     */
    public void initialize(final String loaderJSONFile, final String queryJSONFile,
            final FeatureCodec<T, SOURCE> codec) throws IOException

    {
        mCodec = codec;
        mLoaderJSONFile = loaderJSONFile;
        mQueryJSONFile = queryJSONFile;
        //Read header
        GenomicsDBQueryStream gdbStream = new GenomicsDBQueryStream(loaderJSONFile, queryJSONFile,
                mCodec instanceof BCF2Codec);
        SOURCE source = codec.makeSourceFromStream(gdbStream);
        mFCHeader = codec.readHeader(source);
        //Store sequence names
        mVCFHeader = (VCFHeader)(mFCHeader.getHeaderValue());
        mSequenceNames = new ArrayList<String>(mVCFHeader.getContigLines().size());
        for(final VCFContigHeaderLine contigHeaderLine : mVCFHeader.getContigLines())
            mSequenceNames.add(contigHeaderLine.getID());
        gdbStream.close();
    }

    /**
     * Return the VCF header of the combined gVCF stream
     * @return the VCF header of the combined gVCF stream
     */
    public Object getHeader()
    {
        return mFCHeader.getHeaderValue();
    } 

    /**
     * Return the list of contigs in the combined VCF header
     * @return list of strings of the contig names
     */
    public List<String> getSequenceNames()
    {
        return mSequenceNames;
    }

    public void close() throws IOException
    {
    }

    /**
     * Return an iterator over {@link htsjdk.variant.variantcontext.VariantContext}
     * objects for the specified TileDB array and query configuration
     * @return iterator over {@link htsjdk.variant.variantcontext.VariantContext} objects
     */
    public CloseableTribbleIterator<T> iterator() throws IOException
    {
        return new GenomicsDBFeatureIterator(mLoaderJSONFile, mQueryJSONFile, mCodec);
    }

    /**
     * Return an iterator over {@link htsjdk.variant.variantcontext.VariantContext}
     * objects for the specified TileDB array and queried position
     * @param chr contig name
     * @param start start position (1-based)
     * @param end end position, inclusive (1-based)
     * @return iterator over {@link htsjdk.variant.variantcontext.VariantContext} objects
     */
    public CloseableTribbleIterator<T> query(final String chr, final int start, final int end)
      throws IOException
    {
        return new GenomicsDBFeatureIterator(mLoaderJSONFile, mQueryJSONFile,
          mCodec, chr, start, end);
    }

    /**
     * Iterator over {@link htsjdk.variant.variantcontext.VariantContext} objects.
     * Uses {@link GenomicsDBQueryStream} to obtain combined gVCF records
     * (as BCF2) from TileDB/GenomicsDB
     */
    class GenomicsDBFeatureIterator implements CloseableTribbleIterator<T>
    {
        private FeatureCodec<T, SOURCE> mCodec = null;
        private GenomicsDBQueryStream mStream = null;
        private SOURCE mSource = null;
        private GenomicsDBTimer mTimer = null;
        private boolean mClosedBefore = false;

        /**
         * Constructor
         * @param loaderJSONFile GenomicsDB loader JSON configuration file
         * @param queryJSONFile GenomicsDB query JSON configuration file
         * @param codec FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
         *              and {@link htsjdk.variant.vcf.VCFCodec} are tested
         * @throws IOException when data cannot be read from the stream 
         */
        public GenomicsDBFeatureIterator(final String loaderJSONFile, final String queryJSONFile,
                final FeatureCodec<T, SOURCE> codec) throws IOException
        {
            this(loaderJSONFile, queryJSONFile, codec, "", 0, 0);
        }

        /**
         * Constructor
         * @param loaderJSONFile GenomicsDB loader JSON configuration file
         * @param queryJSONFile GenomicsDB query JSON configuration file
         * @param codec FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
         *              and {@link htsjdk.variant.vcf.VCFCodec} are tested
         * @param chr contig name
         * @param start start position (1-based)
         * @param end end position, inclusive (1-based)
         * @throws IOException when data cannot be read from the stream 
         */
        public GenomicsDBFeatureIterator(final String loaderJSONFile, final String queryJSONFile,
                final FeatureCodec<T, SOURCE> codec,
                final String chr, final int start, final int end) throws IOException
        {
            mCodec = codec;
            boolean readAsBCF = mCodec instanceof BCF2Codec;
            mStream = new GenomicsDBQueryStream(loaderJSONFile, queryJSONFile, chr, start, end,
                    readAsBCF);
            if(readAsBCF) //BCF2 codec provides size of header
                mStream.skip(mFCHeader.getHeaderEnd());
            mSource = codec.makeSourceFromStream(mStream);
            if(!readAsBCF) //VCF Codec must parse out header again since getHeaderEnd() returns 0
                mCodec.readHeader(mSource); //no need to store header anywhere
            mTimer = new GenomicsDBTimer();
        }

        @Override
        public boolean hasNext()
        {
            boolean isDone = (mCodec.isDone(mSource));
            if(isDone)
                close();
            return !isDone;
        }

        @Override
        public T next()
        {
            try
            {
                mTimer.start();
                if(mCodec.isDone(mSource))
                    throw new RuntimeException("No more valid records exist");
                T nextObj = mCodec.decode(mSource);
                mTimer.stop();
                return nextObj;
            }
            catch(IOException e)
            {
                throw new RuntimeException("Unknown codec exception");
            }
        } 

        @Override
        public void close()
        {
            if(!mClosedBefore)
            {
                mTimer.print("GenomicsDB iterator next() timer", System.err);
                mCodec.close(mSource);
                mClosedBefore = true;
            }
        }

        @Override
        public GenomicsDBFeatureIterator iterator()
        {
            return this;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Remove is not supported in Iterators");
        }  
    }
}
