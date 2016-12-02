/**
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

import java.io.IOException;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.lang.reflect.Type;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.ArrayList;
import java.util.List;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ContainerFactory;
import org.json.simple.parser.ParseException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.tribble.CloseableTribbleIterator;
import java.lang.Long;

public final class TestBufferStreamVCF2TileDB
{
    /**
     * Wrapper class to maintain stream state for the test driver program
     * The class can maintains
     * (a) VCFHeader
     * (b) CloseableTribbleIterator<VariantContext>
     * (c) mNextVC the next VariantContext object to be sent to VCF2TileDB iff the buffer interface of VCF2TileDB is used (addBufferStream())
     * and not the Iterator<VariantContext> interface (addSortedVariantContextIterator())
     */
    private static class VCFFileStreamInfo
    {
        public int mStreamIdx = -1;
        public VCFHeader mVCFHeader = null;
        public CloseableTribbleIterator<VariantContext> mIterator = null;
        public VariantContext mNextVC = null;

        /**
         * Constructor
         * @param filename path to VCF file
         */
        public VCFFileStreamInfo(final String fileName) throws IOException
        {
            AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(fileName, new VCFCodec(), false);
            mVCFHeader = (VCFHeader)(reader.getHeader());
            mIterator = reader.iterator();
        }
    }

    /**
     * Factory object to maintain order of keys in simple JSON parsing - use LinkedHashMap
     */
    private static class LinkedHashFactory implements ContainerFactory
    {
        @Override
        public List creatArrayContainer()
        {
            return new ArrayList();
        }

        @Override
        public Map createObjectContainer()
        {
            return new LinkedHashMap();
        }
    }

    /**
     * Sample driver code for testing Java VariantContext write API for GenomicsDB
     * The code shows two ways of using the API (a) Iterator<VariantContext> (b) Directly adding VariantContext objects
     * If "-iterators" is passed as the second argument, method (a) is used.
     */
    public static void main(final String[] args) throws IOException, FileNotFoundException, GenomicsDBException, ParseException
    {
        if(args.length < 2)
        {
            System.err.println("For loading: [-iterators] <loader.json> <stream_name_to_file.json> [bufferCapacity rank lbRowIdx ubRowIdx]");
            System.exit(-1);
        }
        int argsLoaderFileIdx = 0;
        if(args[0].equals("-iterators"))
            argsLoaderFileIdx = 1;
        //Buffer capacity
        long bufferCapacity = (args.length >= argsLoaderFileIdx+3) ? Integer.parseInt(args[argsLoaderFileIdx+2]) : 1024;
        //Specify rank (or partition idx) of this process
        int rank = (args.length >= argsLoaderFileIdx+4) ? Integer.parseInt(args[argsLoaderFileIdx+3]) : 0;
        //Specify smallest row idx from which to start loading - useful for incremental loading into existing array
        long lbRowIdx = (args.length >= argsLoaderFileIdx+5) ? Long.parseLong(args[argsLoaderFileIdx+4]) : 0;
        //Specify largest row idx up to which loading should be performed - for completeness
        long ubRowIdx = (args.length >= argsLoaderFileIdx+6) ? Long.parseLong(args[argsLoaderFileIdx+5]) : Long.MAX_VALUE-1;
        //<loader.json> first arg
        VCF2TileDB loader = new VCF2TileDB(args[argsLoaderFileIdx], rank, lbRowIdx, ubRowIdx);
        //<stream_name_to_file.json> - useful for the driver only
        //JSON file that contains "stream_name": "vcf_file_path" entries
        FileReader mappingReader = new FileReader(args[argsLoaderFileIdx+1]);
        JSONParser parser = new JSONParser();
        LinkedHashMap streamNameToFileName = (LinkedHashMap)parser.parse(mappingReader, new LinkedHashFactory());
        ArrayList<VCFFileStreamInfo> streamInfoVec = new ArrayList<VCFFileStreamInfo>();
        long rowIdx = 0;
        for(Object currObj : streamNameToFileName.entrySet())
        {
            Map.Entry<String, String> entry = (Map.Entry<String, String>)currObj;
            VCFFileStreamInfo currInfo = new VCFFileStreamInfo(entry.getValue());
            //The following 2 lines are not mandatory - use initializeSampleInfoMapFromHeader() iff you know for sure that sample names
            //in the VCF header are globally unique across all streams/files. If not, you have 2 options:
            //(a) specify your own mapping from sample index in the header to SampleInfo object (unique_name, rowIdx) OR
            //(b) specify the mapping in the callset_mapping_file (JSON) and pass null to addSortedVariantContextIterator()
            LinkedHashMap<Integer, VCF2TileDB.SampleInfo> sampleIndexToInfo = new LinkedHashMap<Integer, VCF2TileDB.SampleInfo>();
            rowIdx = VCF2TileDB.initializeSampleInfoMapFromHeader(sampleIndexToInfo, currInfo.mVCFHeader, rowIdx);
            int streamIdx = -1;
            if(args[0].equals("-iterators"))
                streamIdx = loader.addSortedVariantContextIterator(entry.getKey(), currInfo.mVCFHeader, currInfo.mIterator,
                        bufferCapacity, VariantContextWriterBuilder.OutputType.BCF_STREAM, sampleIndexToInfo); //pass sorted VC iterators
            else
                streamIdx = loader.addBufferStream(entry.getKey(), currInfo.mVCFHeader, bufferCapacity,
                        VariantContextWriterBuilder.OutputType.BCF_STREAM, sampleIndexToInfo); //use buffers - VCs will be provided by caller
            currInfo.mStreamIdx = streamIdx;
            streamInfoVec.add(currInfo);
        }
        if(args[0].equals("-iterators"))
        {
            //Much simpler interface if using Iterator<VariantContext>
            loader.importBatch();
            assert loader.isDone();
        }
        else
        {
            //Must be called after all iterators/streams added - no more iterators/streams can be added once
            //this function is called
            loader.setupGenomicsDBImporter();
            //Counts and tracks buffer streams for which new data must be supplied
            //Initialized to all the buffer streams
            int numExhaustedBufferStreams = streamInfoVec.size();
            int[] exhaustedBufferStreamIdxs = new int[numExhaustedBufferStreams];
            for(int i=0;i<numExhaustedBufferStreams;++i)
                exhaustedBufferStreamIdxs[i] = i;
            while(!loader.isDone())
            {
                //Add data for streams that were exhausted in the previous round
                for(int i=0;i<numExhaustedBufferStreams;++i)
                {
                    VCFFileStreamInfo currInfo = streamInfoVec.get(exhaustedBufferStreamIdxs[i]);
                    boolean added = true;
                    while(added && (currInfo.mIterator.hasNext() || currInfo.mNextVC != null))
                    {
                        if(currInfo.mNextVC != null)
                            added = loader.add(currInfo.mNextVC, currInfo.mStreamIdx);
                        if(added)
                            if(currInfo.mIterator.hasNext())
                                currInfo.mNextVC = currInfo.mIterator.next();
                            else
                                currInfo.mNextVC = null;
                    }
                }
                loader.importBatch();
                numExhaustedBufferStreams = (int)loader.getNumExhaustedBufferStreams();
                for(int i=0;i<numExhaustedBufferStreams;++i)
                    exhaustedBufferStreamIdxs[i] = loader.getExhaustedBufferStreamIndex(i);
            }
        }
    }
}
