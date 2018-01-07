/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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
import com.intel.genomicsdb.GenomicsDBCallsetsMapProto;
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
import java.lang.Long;
import gnu.getopt.Getopt;
import gnu.getopt.LongOpt;


public final class TestGenomicsDBImporterWithMergedVCFHeader {

  public enum ArgsIdxEnum
  {
    ARGS_IDX_USE_SAMPLES_IN_ORDER(1000),
    ARGS_IDX_FAIL_IF_UPDATING(1001),
    ARGS_IDX_BATCHSIZE(1002),
    ARGS_IDX_VIDMAP_OUTPUT(1003),
    ARGS_IDX_CALLSET_OUTPUT(1004),
    ARGS_IDX_PASS_AS_BCF(1005),
    ARGS_IDX_VCF_HEADER_OUTPUT(1006),
    ARGS_IDX_AFTER_LAST_ARG_IDX(1007);

    private final int mArgsIdx;
    ArgsIdxEnum(final int idx)
    {
      mArgsIdx = idx;
    }

    int idx()
    {
      return mArgsIdx;
    }
  }
  public static void main(final String[] args)
    throws IOException, GenomicsDBException, ParseException
  {
    final int firstEnumIdx = ArgsIdxEnum.ARGS_IDX_USE_SAMPLES_IN_ORDER.idx();
    LongOpt[] longopts = new LongOpt[10];
    longopts[0] = new LongOpt("use_samples_in_order", LongOpt.NO_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_USE_SAMPLES_IN_ORDER.idx());
    longopts[1] = new LongOpt("fail_if_updating", LongOpt.NO_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_FAIL_IF_UPDATING.idx());
    longopts[2] = new LongOpt("interval", LongOpt.REQUIRED_ARGUMENT, null, 'L');
    longopts[3] = new LongOpt("workspace", LongOpt.REQUIRED_ARGUMENT, null, 'w');
    longopts[4] = new LongOpt("array", LongOpt.REQUIRED_ARGUMENT, null, 'A');
    longopts[5] = new LongOpt("batchsize", LongOpt.REQUIRED_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_BATCHSIZE.idx());
    longopts[6] = new LongOpt("vidmap-output", LongOpt.REQUIRED_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_VIDMAP_OUTPUT.idx());
    longopts[7] = new LongOpt("callset-output", LongOpt.REQUIRED_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_CALLSET_OUTPUT.idx());
    longopts[8] = new LongOpt("pass-as-bcf", LongOpt.NO_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_PASS_AS_BCF.idx());
    longopts[9] = new LongOpt("vcf-header-output", LongOpt.REQUIRED_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_VCF_HEADER_OUTPUT.idx());
    //Arg parsing
    Getopt g = new Getopt("TestGenomicsDBImporterWithMergedVCFHeader", args, "w:A:L:", longopts);
    int c = -1;
    String optarg;
    //Array of enums
    final ArgsIdxEnum[] enumArray = ArgsIdxEnum.values();
    boolean useSamplesInOrder = false;
    boolean failIfUpdating = false;
    String workspace = "";
    String arrayName = "";
    String chromosomeInterval = "";
    int batchSize = 1000000;
    String vidmapOutputFilepath = null;
    String callsetOutputFilepath = null;
    String vcfHeaderOutputFilepath = null;
    boolean passAsVcf = true;
    while ((c = g.getopt()) != -1)
    {
      switch(c)
      {
        case 'w':
          workspace = g.getOptarg();
          break;
        case 'A':
          arrayName = g.getOptarg();
          break;
        case 'L':
          chromosomeInterval = g.getOptarg();
          break;
        default:
          {
            if(c >= firstEnumIdx && c < ArgsIdxEnum.ARGS_IDX_AFTER_LAST_ARG_IDX.idx())
            {
              int offset = c - firstEnumIdx;
              assert offset < enumArray.length;
              switch(enumArray[offset])
              {
                case ARGS_IDX_USE_SAMPLES_IN_ORDER:
                  useSamplesInOrder = true;
                  break;
                case ARGS_IDX_FAIL_IF_UPDATING:
                  failIfUpdating = true;
                  break;
                case ARGS_IDX_BATCHSIZE:
                  batchSize = Integer.parseInt(g.getOptarg());
                  break;
                case ARGS_IDX_VIDMAP_OUTPUT:
                  vidmapOutputFilepath = new String(g.getOptarg());
                  break;
                case ARGS_IDX_CALLSET_OUTPUT:
                  callsetOutputFilepath = new String(g.getOptarg());
                  break;
                case ARGS_IDX_VCF_HEADER_OUTPUT:
                  vcfHeaderOutputFilepath = new String(g.getOptarg());
                  break;
                case ARGS_IDX_PASS_AS_BCF:
                  passAsVcf = false;
                  break;
                default:
                  System.err.println("Unknown command line option "+g.getOptarg()+" - ignored");
                  break;
              }
            }
            else
              System.err.println("Unknown command line option "+g.getOptarg()+" - ignored");
          }
      }
    }
    int numPositionalArgs = args.length - g.getOptind();
    if (numPositionalArgs <= 0
        || arrayName.isEmpty() || workspace.isEmpty()
        || chromosomeInterval.isEmpty()
        ) {
      System.out.println("Usage: ExampleGenomicsDBImporter" + " -L chromosome:interval " +
          "-w genomicsdbworkspace -A arrayname variantfile(s) [--use_samples_in_order --fail_if_updating --batchsize=<N> --vidmap-output <path>]");
      System.exit(-1);
    }

    String[] temp0 = chromosomeInterval.split(":");
    String chromosomeName = temp0[0];
    String[] interval = temp0[1].split("-");
    List<String> files = new ArrayList<>();
    for (int i = g.getOptind(); i < args.length; ++i) {
      files.add(args[i]);
    }

    List<VCFHeader> headers = new ArrayList<>();
    ArrayList<String> sampleNames = new ArrayList<String>();
    Map<String, String> sampleNameToFileName = new LinkedHashMap<String, String>();

    //Get merged header first
    for (String file : files) {
      AbstractFeatureReader<VariantContext, LineIterator> reader =
        AbstractFeatureReader.getFeatureReader(file, new VCFCodec(), false);
      headers.add((VCFHeader) reader.getHeader());
      final String sampleName = ((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0);
      sampleNames.add(sampleName);
      sampleNameToFileName.put(sampleName, file);
      //Hopefully, GC kicks in and frees resources assigned to reader
    }
    Set<VCFHeaderLine> mergedHeader = VCFUtils.smartMergeHeaders(headers, true);

    //This sorts the list sampleNames if !useSamplesInOrder
    //Why you should use this? If you are writing multiple partitions in different machines,
    //you must have consistent ordering of samples across partitions. If file order is different
    //in different processes, then set useSamplesInOrder to false and let the sort in
    //generateSortedCallSetMap ensure consistent ordering across samples
    GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB =
        GenomicsDBImporter.generateSortedCallSetMap(sampleNames, useSamplesInOrder);

    //Write out vidmap if needed
    if(vidmapOutputFilepath != null)
        GenomicsDBImporter.writeVidMapJSONFile(vidmapOutputFilepath, mergedHeader);

    //Write out callset map if needed
    if(callsetOutputFilepath != null)
        GenomicsDBImporter.writeCallsetMapJSONFile(callsetOutputFilepath, callsetMappingPB);

    //write out merged header if needed
    if(vcfHeaderOutputFilepath != null)
        GenomicsDBImporter.writeVcfHeaderFile(vcfHeaderOutputFilepath, mergedHeader);

    //Iterate over sorted sample list in batches
    for(int i=0;i<sampleNames.size();i+=batchSize)
    {
      Map<String, FeatureReader<VariantContext>> map = new LinkedHashMap<>();
      for(int j=i;j<sampleNames.size() && j<i+batchSize;++j)
      {
        final String sampleName = sampleNames.get(j);
        assert sampleNameToFileName.containsKey(sampleName);
        AbstractFeatureReader<VariantContext, LineIterator> reader =
          AbstractFeatureReader.getFeatureReader(sampleNameToFileName.get(sampleName), new VCFCodec(), false);
        assert sampleName.equals(((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0));
        map.put(sampleName, reader);
      }
      GenomicsDBImporter importer = new GenomicsDBImporter(
          map, mergedHeader,
          new ChromosomeInterval(chromosomeName, Integer.parseInt(interval[0]), Integer.parseInt(interval[1])),
          workspace, arrayName, 10000L*sampleNames.size(), 1048576L,
          (long)i, (long)(i+batchSize-1),
          useSamplesInOrder, failIfUpdating,
          0, true, passAsVcf);
      boolean isdone = importer.importBatch();
      assert (isdone);
    }
  }
}
