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

import com.intel.genomicsdb.ChromosomeInterval;
import com.intel.genomicsdb.GenomicsDBException;
import com.intel.genomicsdb.GenomicsDBImporter;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.json.simple.parser.ParseException;

import java.io.IOException;
import java.util.*;

/**
 * Created by kdatta1 on 3/10/17.
 */
public final class TestGenomicsDBImporterWithMergedVCFHeader {

  public static void main(final String[] args)
    throws IOException, GenomicsDBException, ParseException
  {
    if (args.length < 5) {
      System.out.println("Usage: ExampleGenomicsDBImporter" + " chromosome:interval " +
        "genomicsdbworkspace arrayname useSamplesInOrder variantfile(s)");
      System.exit(-1);
    }

    String chromosomeInterval = args[0];
    String[] temp0 = chromosomeInterval.split(":");
    String chromosomeName = temp0[0];
    String[] interval = temp0[1].split("-");
    String workspace = args[1];
    String arrayName = args[2];
    boolean useSamplesInOrder = Boolean.parseBoolean(args[3]);
    List<String> files = new ArrayList<>();
    for (int i = 4; i < args.length; ++i) {
      files.add(args[i]);
    }

    Map<String, FeatureReader<VariantContext>> map = new LinkedHashMap<>();
    List<VCFHeader> headers = new ArrayList<>();

    for (String file : files) {
      AbstractFeatureReader<VariantContext, LineIterator> reader =
        AbstractFeatureReader.getFeatureReader(file, new VCFCodec(), false);
      String sampleName = ((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0);
      headers.add((VCFHeader) reader.getHeader());
      map.put(sampleName, reader);
    }

    Set<VCFHeaderLine> mergedHeader = VCFUtils.smartMergeHeaders(headers, true);

    GenomicsDBImporter importer = new GenomicsDBImporter(
      map, mergedHeader,
      new ChromosomeInterval(chromosomeName, Integer.parseInt(interval[0]), Integer.parseInt(interval[1])),
      workspace, arrayName, 1000L, 1048576L, useSamplesInOrder);
    boolean isdone = importer.importBatch();
    assert (isdone);
  }
}
