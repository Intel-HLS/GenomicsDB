import java.io.IOException;
import java.io.FileReader;
import java.io.FileNotFoundException;
import com.google.gson.Gson;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.ArrayList;
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

import genomicsdb.VCF2TileDB;
import genomicsdb.GenomicsDBException;

public final class TestBufferStreamVCF2TileDB
{
    //Wrapper class to maintain some testing stream state
    private static class VCFFileStreamInfo
    {
        public int mStreamIdx = -1;
        public VCFHeader mVCFHeader = null;
        public CloseableTribbleIterator<VariantContext> mIterator = null;
        public VariantContext mNextVC = null;

        public VCFFileStreamInfo(final String fileName) throws IOException
        {
            AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(fileName, new VCFCodec(), false);
            mVCFHeader = (VCFHeader)(reader.getHeader());
            mIterator = reader.iterator();
        }
    }

    public static void main(final String[] args) throws IOException, FileNotFoundException, GenomicsDBException
    {
        if(args.length < 2)
        {
            System.err.println("For loading: <loader.json> <stream_name_to_file.json> [rank lbRowIdx ubRowIdx]");
            System.exit(-1);
        } 
        //Specify rank (or partition idx) of this process
        int rank = (args.length >= 4) ? Integer.parseInt(args[3]) : 0;
        //Specify smallest row idx from which to start loading - useful for incremental loading into existing array
        long lbRowIdx = (args.length >= 5) ? Long.parseLong(args[4]) : 0;
        //Specify largest row idx up to which loading should be performed - for completeness
        long ubRowIdx = (args.length >= 6) ? Long.parseLong(args[5]) : Long.MAX_VALUE-1;
        //<loader.json> first arg
        VCF2TileDB loader = new VCF2TileDB(args[0], rank, lbRowIdx, ubRowIdx);
        //<stream_name_to_file.json>
        FileReader mappingReader = new FileReader(args[1]);
        Gson gson = new Gson();
        LinkedHashMap<String, String> streamNameToFileName = (LinkedHashMap<String, String>)gson.fromJson(mappingReader, LinkedHashMap.class);
        ArrayList<VCFFileStreamInfo> streamInfoVec = new ArrayList<VCFFileStreamInfo>();
        for(Map.Entry<String, String> entry : streamNameToFileName.entrySet())
        {
            VCFFileStreamInfo currInfo = new VCFFileStreamInfo(entry.getValue());
            int streamIdx = loader.addBufferStream(entry.getKey(), currInfo.mVCFHeader, 10240, VariantContextWriterBuilder.OutputType.BCF_STREAM);
            currInfo.mStreamIdx = streamIdx;
            streamInfoVec.add(currInfo);
        }
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
