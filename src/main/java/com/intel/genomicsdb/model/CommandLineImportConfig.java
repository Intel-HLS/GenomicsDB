package com.intel.genomicsdb.model;

import gnu.getopt.Getopt;
import gnu.getopt.LongOpt;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFUtils;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Consumer;
import java.util.stream.IntStream;

import static java.util.stream.Collectors.toList;

public class CommandLineImportConfig extends ImportConfig {
    //TODO: remove this constant once C++ layer makes use of protobuf structures
    protected static final long DEFAULT_SIZE_PER_COLUMN_PARTITION = 16384L;
    private GenomicsDBImportConfiguration.ImportConfiguration.Builder configurationBuilder =
            GenomicsDBImportConfiguration.ImportConfiguration.newBuilder();
    private List<GenomicsDBImportConfiguration.Partition.Builder> partitionList = new ArrayList<>();
    private String workspace = "";
    private String vcfOutputFilename = "";
    private Consumer<String[]> chromInterToPartitionList = par -> {
        GenomicsDBImportConfiguration.Partition.Builder partitionBuilder = GenomicsDBImportConfiguration.Partition.newBuilder();
        Coordinates.ContigPosition.Builder contigPositionBuilder = Coordinates.ContigPosition.newBuilder();
        Coordinates.GenomicsDBColumn.Builder columnBuilder = Coordinates.GenomicsDBColumn.newBuilder();
        //begin
        contigPositionBuilder.setContig(par[0]).setPosition(Long.parseLong(par[1].split("-")[0]));
        columnBuilder.setContigPosition(contigPositionBuilder.build());
        partitionBuilder.setBegin(columnBuilder.build());
        //end
        contigPositionBuilder.setPosition(Long.parseLong(par[1].split("-")[1]));
        columnBuilder.setContigPosition(contigPositionBuilder.build());
        partitionBuilder.setEnd(columnBuilder.build());
        partitionList.add(partitionBuilder);
    };

    public CommandLineImportConfig(final String command, final String[] commandArgs) {
        Getopt getOpt = new Getopt(command, commandArgs, "w:A:L:", resolveLongOpt());
        //TODO: remove next line once C++ layer makes use of protobuf structures. Making size per column partition explicit.
        configurationBuilder.setSizePerColumnPartition(DEFAULT_SIZE_PER_COLUMN_PARTITION);
        resolveCommandArgs(getOpt);
        partitionList.forEach(partition -> {
                    partition.setWorkspace(this.workspace);
                    if (!vcfOutputFilename.isEmpty()) partition.setVcfOutputFilename(vcfOutputFilename);
                });
        partitionList.forEach(partition -> configurationBuilder.addColumnPartitions(partition));
        this.setImportConfiguration(configurationBuilder.build());
        this.validateChromosomeIntervals();
        int numPositionalArgs = commandArgs.length - getOpt.getOptind();
        if (numPositionalArgs <= 0
                || this.getImportConfiguration().getColumnPartitions(0).getWorkspace().isEmpty()
                || this.getImportConfiguration().getColumnPartitionsList().isEmpty()) {
            throwIllegalArgumentException();
        }
        List<String> files = IntStream.range(getOpt.getOptind(), commandArgs.length).mapToObj(
                i -> commandArgs[i]).collect(toList());
        this.resolveHeaders(files);
        this.setSampleToReaderMapCreator(this::createSampleToReaderMap);
    }

    private LongOpt[] resolveLongOpt() {
        LongOpt[] longopts = new LongOpt[11];
        longopts[0] = new LongOpt("use_samples_in_order", LongOpt.NO_ARGUMENT, null,
                ArgsIdxEnum.ARGS_IDX_USE_SAMPLES_IN_ORDER.idx());
        longopts[1] = new LongOpt("fail_if_updating", LongOpt.NO_ARGUMENT, null,
                ArgsIdxEnum.ARGS_IDX_FAIL_IF_UPDATING.idx());
        longopts[2] = new LongOpt("interval", LongOpt.REQUIRED_ARGUMENT, null, 'L');
        longopts[3] = new LongOpt("workspace", LongOpt.REQUIRED_ARGUMENT, null, 'w');
        longopts[4] = new LongOpt("batchsize", LongOpt.REQUIRED_ARGUMENT, null,
                ArgsIdxEnum.ARGS_IDX_BATCHSIZE.idx());
        longopts[5] = new LongOpt("vidmap-output", LongOpt.REQUIRED_ARGUMENT, null,
                ArgsIdxEnum.ARGS_IDX_VIDMAP_OUTPUT.idx());
        longopts[6] = new LongOpt("callset-output", LongOpt.REQUIRED_ARGUMENT, null,
                ArgsIdxEnum.ARGS_IDX_CALLSET_OUTPUT.idx());
        longopts[7] = new LongOpt("pass-as-bcf", LongOpt.NO_ARGUMENT, null,
                ArgsIdxEnum.ARGS_IDX_PASS_AS_BCF.idx());
        longopts[8] = new LongOpt("vcf-header-output", LongOpt.REQUIRED_ARGUMENT, null,
                ArgsIdxEnum.ARGS_IDX_VCF_HEADER_OUTPUT.idx());
        longopts[9] = new LongOpt("size_per_column_partition", LongOpt.REQUIRED_ARGUMENT, null,
                ArgsIdxEnum.ARGS_IDX_SIZE_PER_COLUMN_PARTITION.idx());
        longopts[10] = new LongOpt("segment_size", LongOpt.REQUIRED_ARGUMENT, null,
                ArgsIdxEnum.ARGS_IDX_SEGMENT_SIZE.idx());
        return longopts;
    }

    private void resolveCommandArgs(final Getopt commandArgs) {
        int c;
        final int firstEnumIdx = ArgsIdxEnum.ARGS_IDX_USE_SAMPLES_IN_ORDER.idx();
        final ArgsIdxEnum[] enumArray = ArgsIdxEnum.values();
        while ((c = commandArgs.getopt()) != -1) {
            switch (c) {
                case 'w':
                    workspace = commandArgs.getOptarg();
                    break;
                case 'L':
                    chromInterToPartitionList.accept(commandArgs.getOptarg().split(":"));
                    break;
                default: {
                    if (c >= firstEnumIdx && c < ArgsIdxEnum.ARGS_IDX_AFTER_LAST_ARG_IDX.idx()) {
                        int offset = c - firstEnumIdx;
                        assert offset < enumArray.length;
                        switch (enumArray[offset]) {
                            case ARGS_IDX_USE_SAMPLES_IN_ORDER:
                                this.setUseSamplesInOrder(true);
                                break;
                            case ARGS_IDX_FAIL_IF_UPDATING:
                                configurationBuilder.setFailIfUpdating(true);
                                break;
                            case ARGS_IDX_BATCHSIZE:
                                setBatchSize(Integer.parseInt(commandArgs.getOptarg()));
                                break;
                            case ARGS_IDX_VIDMAP_OUTPUT:
                                this.setOutputVidmapJsonFile(commandArgs.getOptarg());
                                break;
                            case ARGS_IDX_CALLSET_OUTPUT:
                                this.setOutputCallsetmapJsonFile(commandArgs.getOptarg());
                                break;
                            case ARGS_IDX_VCF_HEADER_OUTPUT:
                                vcfOutputFilename = commandArgs.getOptarg();
                                break;
                            case ARGS_IDX_PASS_AS_BCF:
                                setPassAsVcf(false);
                                break;
                            case ARGS_IDX_SIZE_PER_COLUMN_PARTITION:
                                configurationBuilder.setSizePerColumnPartition(Long.parseLong(commandArgs.getOptarg()));
                                break;
                            case ARGS_IDX_SEGMENT_SIZE:
                                configurationBuilder.setSegmentSize(Long.parseLong(commandArgs.getOptarg()));
                                break;
                            default:
                                throwIllegalArgumentException(commandArgs);
                                return;
                        }
                    } else {
                        throwIllegalArgumentException(commandArgs);
                    }
                }
            }
        }
    }

    private void throwIllegalArgumentException() {
        throw new IllegalArgumentException("Invalid usage. Correct way of using arguments: -L chromosome:interval " +
                "-w genomicsdbworkspace --size_per_column_partition 10000 --segment_size 1048576 variantfile(s) " +
                "[--use_samples_in_order --fail_if_updating --batchsize=<N> --vidmap-output <path>]");
    }

    private void throwIllegalArgumentException(Getopt commandArgs) {
        throw new IllegalArgumentException("Unknown command line option " +
                commandArgs.getOptarg() + " - ignored");
    }

    private void resolveHeaders(final List<String> files) {
        List<VCFHeader> headers = new ArrayList<>();
        Map<String, Path> sampleNameToVcfPath = new LinkedHashMap<>();

        //Get merged header first
        for (String file : files) {
            AbstractFeatureReader<VariantContext, LineIterator> reader =
                    AbstractFeatureReader.getFeatureReader(file, new VCFCodec(), false);
            headers.add((VCFHeader) reader.getHeader());
            final String sampleName = ((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0);
            sampleNameToVcfPath.put(sampleName, Paths.get(file));
            //Hopefully, GC kicks in and frees resources assigned to reader
        }

        this.setMergedHeader(VCFUtils.smartMergeHeaders(headers, true));
        this.setSampleNameToVcfPath(sampleNameToVcfPath);
    }

    private Map<String, FeatureReader<VariantContext>> createSampleToReaderMap(
            final Map<String, Path> sampleNameToVcfPath, final int batchSize, final int index) {
        final Map<String, FeatureReader<VariantContext>> sampleToReaderMap = new LinkedHashMap<>();
        for (int j = index; j < sampleNameToVcfPath.size() && j < index + batchSize; ++j) {
            final String sampleName = sampleNameToVcfPath.keySet().toArray()[j].toString(); //TODO: fix this.
            assert sampleNameToVcfPath.containsKey(sampleName);
            AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                    sampleNameToVcfPath.get(sampleName).toAbsolutePath().toString(), new VCFCodec(), false);
            assert sampleName.equals(((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0));
            sampleToReaderMap.put(sampleName, reader);
        }
        return sampleToReaderMap;
    }

    public enum ArgsIdxEnum {
        ARGS_IDX_USE_SAMPLES_IN_ORDER(1000),
        ARGS_IDX_FAIL_IF_UPDATING(1001),
        ARGS_IDX_BATCHSIZE(1002),
        ARGS_IDX_VIDMAP_OUTPUT(1003),
        ARGS_IDX_CALLSET_OUTPUT(1004),
        ARGS_IDX_PASS_AS_BCF(1005),
        ARGS_IDX_VCF_HEADER_OUTPUT(1006),
        ARGS_IDX_SIZE_PER_COLUMN_PARTITION(1007),
        ARGS_IDX_SEGMENT_SIZE(1008),
        ARGS_IDX_AFTER_LAST_ARG_IDX(1009);
        private final int mArgsIdx;

        ArgsIdxEnum(final int idx) {
            mArgsIdx = idx;
        }

        public int idx() {
            return mArgsIdx;
        }
    }
}
