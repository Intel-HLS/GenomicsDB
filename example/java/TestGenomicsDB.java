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
import com.intel.genomicsdb.GenomicsDBFeatureReader;
import com.intel.genomicsdb.GenomicsDBImporter;

public final class TestGenomicsDB
{
  public static void main(final String[] args) throws IOException
  {
    if(args.length < 2)
    {
      System.err.println("Usage:\n\tFor querying: -query <loader.json> [<query.json> |"
        +"<workspace> <array> <reference_genome> <template_VCF_header>"
        +" [<chr> <start> <end>] ]\n"
        +"\tFor loading: -load <loader.json> [rank lbRowIdx ubRowIdx]");
      System.exit(-1);
    }
    //Do query
    if(args[0].equals("-query"))
    {
      GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> reader = null;
      if(args.length == 3)
        //query with <loader.json> <query.json> - required positions may be
        //specified in <query.json>
        //If no positions are specified in <query.json>, the whole array will be scanned
        reader = new GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream>
          (args[1], args[2], new BCF2Codec());
      else
      if(args.length >= 6)
        //query with <loader.json> <workspace> <array> <reference_genome>
        //<template_VCF_header>
        reader = new GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream>
          (args[1], args[2], args[3], args[4], args[5], new BCF2Codec());
      final VariantContextWriter writer =
        new VariantContextWriterBuilder().setOutputVCFStream(System.out).unsetOption(
          Options.INDEX_ON_THE_FLY).build();
      writer.writeHeader((VCFHeader)(reader.getHeader()));
      if(args.length == 6 || args.length == 3) //chr not specified on command line
      {
        CloseableTribbleIterator<VariantContext> gdbIterator = reader.iterator();
        while(gdbIterator.hasNext())
          writer.add(gdbIterator.next());
        gdbIterator.close();
      }
      else
      if(args.length >= 9)
        //chr,start,end specified on the command line - scan only the required positions
        for(final VariantContext record : reader.query(args[6],
          Integer.parseInt(args[7]), Integer.parseInt(args[8])))
          writer.add(record);
    }
    else
    if(args[0].equals("-load")) //load data
    {
      //<loader.json>
      GenomicsDBImporter loader = new GenomicsDBImporter(args[1]);
      //Specify rank (or partition idx) of this process
      int rank = (args.length >= 3) ? Integer.parseInt(args[2]) : 0;
      //Specify smallest row idx from which to start loading - useful for
      //incremental loading into existing array
      long lbRowIdx = (args.length >= 4) ? Long.parseLong(args[3]) : 0;
      //Specify largest row idx up to which loading should be performed - for completeness
      long ubRowIdx = (args.length >= 5) ? Long.parseLong(args[4]) : Long.MAX_VALUE-1;
      loader.write(rank, lbRowIdx, ubRowIdx);
    }
  }
}
