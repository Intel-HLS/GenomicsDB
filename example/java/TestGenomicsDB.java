import java.io.IOException;

import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.vcf.VCFHeader;

import genomicsdb.GenomicsDBFeatureReader;

public final class TestGenomicsDB
{
    public static void main(final String[] args) throws IOException
    {
        if(args.length < 2)
        {
            System.err.println("Usage: <loader.json> [<query.json> | <workspace> <array> <reference_genome> <template_VCF_header>]");
            System.exit(-1);
        }
        GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> reader = null;
        if(args.length == 2)
            reader = new GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream>(args[0], args[1], new BCF2Codec());
        else
            if(args.length >= 5)
                reader = new GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream>(args[0], args[1], args[2], args[3], args[4],
                        new BCF2Codec());
        final VariantContextWriter writer = new VariantContextWriterBuilder().setOutputVCFStream(System.out).unsetOption(Options.INDEX_ON_THE_FLY).build();
        writer.writeHeader((VCFHeader)(reader.getHeader()));
        if(args.length == 5 || args.length == 2)
            for(final VariantContext record : reader.iterator())
                writer.add(record);
        else
            if(args.length >= 8)
                for(final VariantContext record : reader.query(args[5], Integer.parseInt(args[6]), Integer.parseInt(args[7])))
                    writer.add(record);
    }
}
