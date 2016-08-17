import java.io.IOException;

import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.tribble.CloseableTribbleIterator;
import java.lang.Long;

import genomicsdb.GenomicsDBFeatureReader;
import genomicsdb.VCF2TileDB;

public final class TestGenomicsDB
{
    public static void main(final String[] args) throws IOException
    {
        if(args.length < 2)
        {
            System.err.println("Usage:\n\tFor querying: -query <loader.json> [<query.json> | <workspace> <array> <reference_genome> <template_VCF_header>]\n"+
                    "\tFor loading: -load <loader.json> [rank lbSampleIdx ubSampleIdx]");
            System.exit(-1);
        }
        if(args[0].equals("-query"))
        {
            GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> reader = null;
            if(args.length == 3)
                reader = new GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream>(args[1], args[2], new BCF2Codec());
            else
                if(args.length >= 6)
                    reader = new GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream>(args[1], args[2], args[3], args[4], args[5],
                            new BCF2Codec());
            final VariantContextWriter writer = new VariantContextWriterBuilder().setOutputVCFStream(System.out).unsetOption(Options.INDEX_ON_THE_FLY).build();
            writer.writeHeader((VCFHeader)(reader.getHeader()));
            if(args.length == 6 || args.length == 3)
            {
                CloseableTribbleIterator<VariantContext> gdbIterator = reader.iterator();
                while(gdbIterator.hasNext())
                    writer.add(gdbIterator.next());
                gdbIterator.close();
            }
            else
                if(args.length >= 9)
                    for(final VariantContext record : reader.query(args[6], Integer.parseInt(args[7]), Integer.parseInt(args[8])))
                        writer.add(record);
        }
        else
            if(args[0].equals("-load"))
            {
                VCF2TileDB loader = new VCF2TileDB(args[1]);
                int rank = (args.length >= 3) ? Integer.parseInt(args[2]) : 0;
                long lbSampleIdx = (args.length >= 4) ? Long.parseLong(args[3]) : 0;
                long ubSampleIdx = (args.length >= 5) ? Long.parseLong(args[4]) : Long.MAX_VALUE-1;
                loader.write(rank, lbSampleIdx, ubSampleIdx);
            }
    }
}
