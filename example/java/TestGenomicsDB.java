/*
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

import java.io.IOException;

import com.googlecode.protobuf.format.JsonFormat;
import com.intel.genomicsdb.GenomicsDBExportConfiguration;
import com.intel.genomicsdb.GenomicsDBImporter;
import com.intel.genomicsdb.reader.GenomicsDBFeatureReader;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.tribble.CloseableTribbleIterator;
import java.lang.Long;
import java.lang.Integer;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Optional;

import gnu.getopt.Getopt;
import gnu.getopt.LongOpt;

public final class TestGenomicsDB
{
  private static String readFile(String path, Charset encoding) throws IOException {
    byte[] encoded = Files.readAllBytes(Paths.get(path));
    return new String(encoded, encoding);
  }

  private static GenomicsDBExportConfiguration.ExportConfiguration resolveExportConfigurationFromJsonString(String jsonString) throws JsonFormat.ParseException {
    GenomicsDBExportConfiguration.ExportConfiguration.Builder exportConfigurationBuilder =
            GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();
    JsonFormat.merge(jsonString, exportConfigurationBuilder);
    return exportConfigurationBuilder.build();
  }

  public static <SourceType, CodecType extends FeatureCodec<VariantContext, SourceType>> void runQuery(
          CodecType codec,
          String[] args, int optind, int numPositionalArgs,
          final String loaderJSONFile,
          final String workspace, final String array, final String referenceGenome, final String templateVCFHeader,
          final String chromosome, final int chrBegin, final int chrEnd,
          final boolean countOnly) throws IOException
  {
    GenomicsDBFeatureReader<VariantContext, SourceType> reader;
    GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration;
    if(numPositionalArgs == 1)
    {
      //query with [-l <loader.json>] <query.json> - required positions may be
      //specified in <query.json>
      //Vid and callset mapping file may be specified in the query JSON
      //If no positions are specified in <query.json>, the whole array will be scanned
      String queryJSONFilePath = args[optind];
      String queryJsonFileContent = readFile(queryJSONFilePath, Charset.defaultCharset());
      exportConfiguration = resolveExportConfigurationFromJsonString(queryJsonFileContent);
      reader = new GenomicsDBFeatureReader<>(exportConfiguration, codec, Optional.of(loaderJSONFile));
    }
    else
    {
      //Must specify workspace, array, reference_genome, loader JSON
      //query with <workspace> <array> <reference_genome>
      //<template_VCF_header>
      assert(!loaderJSONFile.isEmpty());
      assert(!workspace.isEmpty());
      assert(!array.isEmpty());
      assert(!referenceGenome.isEmpty());
      exportConfiguration = GenomicsDBExportConfiguration.ExportConfiguration.newBuilder()
              .setWorkspace(workspace)
              .setArray(array)
              .setReferenceGenome(referenceGenome)
              .setVcfHeaderFilename(templateVCFHeader).build();
      reader = new GenomicsDBFeatureReader<>(exportConfiguration, codec, Optional.of(loaderJSONFile));
    }
    final VariantContextWriter writer =
            new VariantContextWriterBuilder().setOutputVCFStream(System.out).unsetOption(
                    Options.INDEX_ON_THE_FLY).build();
    if(!countOnly)
      writer.writeHeader((VCFHeader)(reader.getHeader()));
    if(chromosome.isEmpty()) //chr not specified on command line
    {
      CloseableTribbleIterator<VariantContext> gdbIterator = reader.iterator();
      if(countOnly)
      {
        long counter = 0;
        while(gdbIterator.hasNext())
        {
          ++counter;
          gdbIterator.next();
        }
        System.err.println("#VC objects "+counter);
      }
      else
      {
        while(gdbIterator.hasNext())
          writer.add(gdbIterator.next());
      }
      gdbIterator.close();
    }
    else
    {
      //chr,start,end specified on the command line - scan only the required positions
      CloseableTribbleIterator<VariantContext> gdbIterator = reader.query(chromosome,
              chrBegin, chrEnd);
      if(countOnly)
      {
        long counter = 0;
        while(gdbIterator.hasNext())
        {
          ++counter;
          gdbIterator.next();
        }
        System.err.println("#VC objects "+counter);
      }
      else
      {
        while(gdbIterator.hasNext())
          writer.add(gdbIterator.next());
      }
      gdbIterator.close();
    }
  }

  public enum ArgsIdxEnum
  {
    ARGS_IDX_DO_QUERY(1000),
    ARGS_IDX_DO_LOAD(1001),
    ARGS_IDX_REFERENCE_GENOME(1002),
    ARGS_IDX_TEMPLATE_VCF_HEADER(1003),
    ARGS_IDX_LB_ROW_IDX(1004),
    ARGS_IDX_UB_ROW_IDX(1005),
    ARGS_IDX_CHROMOSOME(1006),
    ARGS_IDX_BEGIN(1007),
    ARGS_IDX_END(1008),
    ARGS_IDX_COUNT_ONLY(1009),
    ARGS_IDX_PASS_AS_VCF(1010),
    ARGS_IDX_AFTER_LAST_ARG_IDX(1011);

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

  public static void main(final String[] args) throws IOException
  {
    int firstEnumIdx = ArgsIdxEnum.ARGS_IDX_DO_QUERY.idx();
    LongOpt[] longopts = new LongOpt[15];
    longopts[0] = new LongOpt("query", LongOpt.NO_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_DO_QUERY.idx());
    longopts[1] = new LongOpt("load", LongOpt.NO_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_DO_LOAD.idx());
    //Specify rank (or partition idx) of this process
    longopts[2] = new LongOpt("rank", LongOpt.REQUIRED_ARGUMENT, null, 'r');
    longopts[5] = new LongOpt("workspace", LongOpt.REQUIRED_ARGUMENT, null, 'w');
    longopts[6] = new LongOpt("array", LongOpt.REQUIRED_ARGUMENT, null, 'A');
    longopts[3] = new LongOpt("reference_genome", LongOpt.REQUIRED_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_REFERENCE_GENOME.idx());
    longopts[4] = new LongOpt("template_vcf_header", LongOpt.REQUIRED_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_TEMPLATE_VCF_HEADER.idx());
    //Specify smallest row idx from which to start loading - useful for
    //incremental loading into existing array
    longopts[7] = new LongOpt("lb_row_idx", LongOpt.REQUIRED_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_LB_ROW_IDX.idx());
    //Specify largest row idx up to which loading should be performed - for completeness
    longopts[8] = new LongOpt("ub_row_idx", LongOpt.REQUIRED_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_UB_ROW_IDX.idx());
    longopts[9] = new LongOpt("chromosome", LongOpt.REQUIRED_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_CHROMOSOME.idx());
    longopts[10] = new LongOpt("begin", LongOpt.REQUIRED_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_BEGIN.idx());
    longopts[11] = new LongOpt("end", LongOpt.REQUIRED_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_END.idx());
    longopts[12] = new LongOpt("loader_json_file", LongOpt.REQUIRED_ARGUMENT, null, 'l');
    longopts[13] = new LongOpt("count_only", LongOpt.NO_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_COUNT_ONLY.idx());
    longopts[14] = new LongOpt("pass_as_vcf", LongOpt.NO_ARGUMENT, null, ArgsIdxEnum.ARGS_IDX_PASS_AS_VCF.idx());
    if(args.length < 2)
    {
      System.err.println("Usage:\n\tFor querying: --query <loader.json> [<query.json> |"
              +"--workspace=<workspace> --array=<array> --reference_genome=<reference_genome> [--template_vcf_header=<template_VCF_header>]"
              +" [--chromosome=<chr> --begin=<start> --end=<end>] ]\n"
              +"\tFor loading: --load <loader.json> [--rank=rank --lb_row_idx=lbRowIdx --ub_row_idx=ubRowIdx]");
      System.exit(-1);
    }
    boolean doQuery = false;
    boolean doLoad = false;
    String loaderJSONFile = "";
    String queryJSONFile = "";
    int rank = 0;
    String workspace = "";
    String array = "";
    String referenceGenome = "";
    String templateVCFHeader = "";
    long lbRowIdx = 0;
    long ubRowIdx = Long.MAX_VALUE-1;
    String chromosome = "";
    int chrBegin = 1;
    int chrEnd = Integer.MAX_VALUE-1;
    boolean countOnly = false;
    boolean passAsVCF = false;
    //Arg parsing
    Getopt g = new Getopt("TestGenomicsDB", args, "w:A:r:l:", longopts);
    int c = -1;
    String optarg;
    //Array of enums
    final ArgsIdxEnum[] enumArray = ArgsIdxEnum.values();
    while ((c = g.getopt()) != -1)
    {
      switch(c)
      {
        case 'r':
          rank = Integer.parseInt(g.getOptarg());
          break;
        case 'w':
          workspace = g.getOptarg();
          break;
        case 'A':
          array = g.getOptarg();
          break;
        case 'l':
          loaderJSONFile = g.getOptarg();
          break;
        default:
        {
          if(c >= firstEnumIdx && c < ArgsIdxEnum.ARGS_IDX_AFTER_LAST_ARG_IDX.idx())
          {
            int offset = c - firstEnumIdx;
            assert offset < enumArray.length;
            switch(enumArray[offset])
            {
              case ARGS_IDX_DO_QUERY:
                doQuery = true;
                break;
              case ARGS_IDX_DO_LOAD:
                doLoad = true;
                break;
              case ARGS_IDX_REFERENCE_GENOME:
                referenceGenome = g.getOptarg();
                break;
              case ARGS_IDX_TEMPLATE_VCF_HEADER:
                templateVCFHeader = g.getOptarg();
                break;
              case ARGS_IDX_LB_ROW_IDX:
                lbRowIdx = Long.parseLong(g.getOptarg());
                break;
              case ARGS_IDX_UB_ROW_IDX:
                ubRowIdx = Long.parseLong(g.getOptarg());
                break;
              case ARGS_IDX_CHROMOSOME:
                chromosome = g.getOptarg();
                break;
              case ARGS_IDX_BEGIN:
                chrBegin = Integer.parseInt(g.getOptarg());
                break;
              case ARGS_IDX_END:
                chrEnd = Integer.parseInt(g.getOptarg());
                break;
              case ARGS_IDX_COUNT_ONLY:
                countOnly = true;
                break;
              case ARGS_IDX_PASS_AS_VCF:
                passAsVCF = true;
                break;
              default:
                System.err.println("Unknown command line option "+g.getOptarg()+" - ignored");
                break;
            }
          }
          else
            System.err.println("Unknown command line option "+g.getOptarg()+" - ignored");
          break;
        }
      }
    }
    //Do either query or load but not both
    assert (!doQuery || !doLoad);
    if(!doLoad && !doQuery)
      doQuery = true;
    int numPositionalArgs = args.length - g.getOptind();
    //Do query
    if(doQuery)
    {
      if(passAsVCF)
        TestGenomicsDB.runQuery(
                new VCFCodec(),
                args, g.getOptind(), numPositionalArgs,
                loaderJSONFile,
                workspace, array, referenceGenome, templateVCFHeader,
                chromosome, chrBegin, chrEnd,
                countOnly);
      else
        TestGenomicsDB.runQuery(
                new BCF2Codec(),
                args, g.getOptind(), numPositionalArgs,
                loaderJSONFile,
                workspace, array, referenceGenome, templateVCFHeader,
                chromosome, chrBegin, chrEnd,
                countOnly);
    }
    else
    {
      if(numPositionalArgs == 1 && loaderJSONFile.isEmpty())
        loaderJSONFile = args[g.getOptind()];
      //<loader.json>
      GenomicsDBImporter loader = new GenomicsDBImporter(loaderJSONFile);
      loader.write(rank, lbRowIdx, ubRowIdx);
    }
  }
}