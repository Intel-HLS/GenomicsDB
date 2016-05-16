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

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFContigHeaderLine;

import java.io.InputStream;
import java.io.IOException;
import java.io.File;
import java.io.ByteArrayInputStream;
import java.io.FileWriter;
import java.util.List;
import java.util.ArrayList;
import java.lang.RuntimeException;

/**
 * A reader for GenomicsDB that provides an iterator over VariantContext objects
 */

public class GenomicsDBFeatureReader<T extends Feature, SOURCE> implements FeatureReader<T>
{
    private String mLoaderJSONFile = null;
    private String mQueryJSONFile = null;
    private FeatureCodec<T, SOURCE> mCodec = null;
    private SOURCE mSource;
    protected FeatureCodecHeader mFCHeader;
    private VCFHeader mVCFHeader = null;
    private ArrayList<String> mSequenceNames;

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
        queryJSON += indentString + "\"vcf_header_filename\": \""+templateVCFHeaderFilename+"\"\n";
        queryJSON += "}\n";
        File tmpQueryJSONFile = File.createTempFile("queryJSON", ".json");
        tmpQueryJSONFile.deleteOnExit();
        FileWriter fptr = new FileWriter(tmpQueryJSONFile);
        fptr.write(queryJSON);
        fptr.close();
        initialize(loaderJSONFile, tmpQueryJSONFile.getAbsolutePath(), codec);
    }

    public GenomicsDBFeatureReader(final String loaderJSONFile, final String queryJSONFile,
            final FeatureCodec<T, SOURCE> codec) throws IOException
    {
        initialize(loaderJSONFile, queryJSONFile, codec);
    }

    public void initialize(final String loaderJSONFile, final String queryJSONFile,
            final FeatureCodec<T, SOURCE> codec) throws IOException

    {
        mCodec = codec;
        mLoaderJSONFile = loaderJSONFile;
        mQueryJSONFile = queryJSONFile;
        //Read header
        GenomicsDBQueryOutputStream gdbStream = new GenomicsDBQueryOutputStream(loaderJSONFile, queryJSONFile);
        mSource = codec.makeSourceFromStream(gdbStream);
        mFCHeader = codec.readHeader(mSource);
        //Store sequence names
        mVCFHeader = (VCFHeader)(mFCHeader.getHeaderValue());
        mSequenceNames = new ArrayList<String>(mVCFHeader.getContigLines().size());
        for(final VCFContigHeaderLine contigHeaderLine : mVCFHeader.getContigLines())
            mSequenceNames.add(contigHeaderLine.getID());
        gdbStream.close();
    }

    public Object getHeader()
    {
        return mFCHeader.getHeaderValue();
    } 

    public List<String> getSequenceNames()
    {
        return mSequenceNames;
    }

    public void close() throws IOException
    {
    }

    public CloseableTribbleIterator<T> iterator() throws IOException
    {
        return new GenomicsDBFeatureIterator(mLoaderJSONFile, mQueryJSONFile, mCodec);
    }

    public CloseableTribbleIterator<T> query(final String chr, final int start, final int end) throws IOException
    {
        return new GenomicsDBFeatureIterator(mLoaderJSONFile, mQueryJSONFile, mCodec, chr, start, end);
    }

    class GenomicsDBFeatureIterator implements CloseableTribbleIterator<T>
    {
        private FeatureCodec<T, SOURCE> mCodec = null;
        private GenomicsDBQueryOutputStream mStream = null;
        private SOURCE mSource = null;

        public GenomicsDBFeatureIterator(final String loaderJSONFile, final String queryJSONFile,
                final FeatureCodec<T, SOURCE> codec) throws IOException
        {
            this(loaderJSONFile, queryJSONFile, codec, "", 0, 0);
        }

        public GenomicsDBFeatureIterator(final String loaderJSONFile, final String queryJSONFile,
                final FeatureCodec<T, SOURCE> codec,
                final String chr, final int start, final int end) throws IOException
        {
            mCodec = codec;
            mStream = new GenomicsDBQueryOutputStream(loaderJSONFile, queryJSONFile, chr, start, end);
            mStream.skip(mFCHeader.getHeaderEnd());
            mSource = codec.makeSourceFromStream(mStream);
        }

        public boolean hasNext()
        {
            return !(mCodec.isDone(mSource));
        }

        public T next()
        {
            try
            {
                if(mCodec.isDone(mSource))
                    throw new RuntimeException("No more valid records exist");
                return mCodec.decode(mSource);
            }
            catch(IOException e)
            {
                throw new RuntimeException("Unknown codec exception");
            }
        } 

        public void close()
        {
            mCodec.close(mSource);
        }

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
