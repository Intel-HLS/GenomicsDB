/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2018 Intel Corporation
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

package com.intel.genomicsdb.model;

import com.intel.genomicsdb.ChromosomeInterval;
import gnu.getopt.Getopt;
import gnu.getopt.LongOpt;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.IntStream;

import static java.util.stream.Collectors.toList;

public class BaseImportConfig {
    //TODO: too many parameters, need to be re-grouped or reduced based on protobuf classes
    private List<ChromosomeInterval> chromosomeIntervalList = new ArrayList<>();
    private String workspace;
    private Long vcfBufferSizePerColumnPartition;
    private Long segmentSize;
    private Long lbRowIdx;
    private Long ubRowIdx;
    private boolean useSamplesInOrderProvided = false;
    private boolean failIfUpdating = false;
    private int rank;
    private boolean validateSampleToReaderMap;
    private boolean passAsVcf = true;
    // Extra params
    private int batchSize = 1000000;
    private String vidmapOutputFilepath;
    private String callsetOutputFilepath;
    private String vcfHeaderOutputFilepath;
    private List<String> files;

    public BaseImportConfig(final List<ChromosomeInterval> chromosomeIntervalList,
                            final String workspace,
                            final Long vcfBufferSizePerColumnPartition,
                            final Long segmentSize,
                            final Long lbRowIdx,
                            final Long ubRowIdx,
                            final boolean useSamplesInOrderProvided,
                            final boolean failIfUpdating,
                            final int rank,
                            final boolean validateSampleToReaderMap,
                            final boolean passAsVcf,
                            final int batchSize,
                            final String vidmapOutputFilepath,
                            final String callsetOutputFilepath,
                            final String vcfHeaderOutputFilepath,
                            final List<String> files) {
        this.validateChromosomeIntervals(chromosomeIntervalList);
        this.validate(workspace, (value) -> value != null && !value.isEmpty(), "Workspace name cannot be empty.");
        this.setChromosomeIntervalList(chromosomeIntervalList);
        this.setWorkspace(workspace);
        this.setVcfBufferSizePerColumnPartition(vcfBufferSizePerColumnPartition);
        this.setSegmentSize(segmentSize);
        this.setLbRowIdx(lbRowIdx);
        this.setUbRowIdx(ubRowIdx);
        this.setUseSamplesInOrderProvided(useSamplesInOrderProvided);
        this.setFailIfUpdating(failIfUpdating);
        this.setRank(rank);
        this.setValidateSampleToReaderMap(validateSampleToReaderMap);
        this.setPassAsVcf(passAsVcf);
        this.setBatchSize(batchSize);
        this.setVidmapOutputFilepath(vidmapOutputFilepath);
        this.setCallsetOutputFilepath(callsetOutputFilepath);
        this.setVcfHeaderOutputFilepath(vcfHeaderOutputFilepath);
        this.setFiles(files);
    }

    private <T> void validate(T attribute, Function<T, Boolean> validation, String errorMessage) {
        if(!validation.apply(attribute)) throw new IllegalArgumentException(errorMessage);
    }

    private void validateChromosomeIntervals(List<ChromosomeInterval> chromosomeIntervalList) {
        if (isThereChromosomeIntervalIntersection(chromosomeIntervalList))
            throw new IllegalArgumentException("There are multiple intervals sharing same value. This is not allowed. " +
                    "Intervals should be defined without intersections.");
    }

    public BaseImportConfig(final String command, final String[] commandArgs) {
        Getopt getOpt = new Getopt(command, commandArgs, "w:A:L:", resolveLongOpt());
        resolveCommandArgs(getOpt);
        this.validateChromosomeIntervals(this.chromosomeIntervalList);
        int numPositionalArgs = commandArgs.length - getOpt.getOptind();
        if (numPositionalArgs <= 0 || this.getWorkspace().isEmpty() || this.getChromosomeIntervalList().isEmpty()) {
            throw new IllegalArgumentException("Invalid usage. Correct way of using arguments: -L chromosome:interval " +
                    "-w genomicsdbworkspace variantfile(s) [--use_samples_in_order --fail_if_updating " +
                    "--batchsize=<N> --vidmap-output <path>]");
        }

        setFiles(IntStream.range(getOpt.getOptind(), commandArgs.length).mapToObj(i -> commandArgs[i]).collect(toList()));
    }

    private LongOpt[] resolveLongOpt() {
        LongOpt[] longopts = new LongOpt[10];
        longopts[0] = new LongOpt("use_samples_in_order", LongOpt.NO_ARGUMENT, null, BaseImportConfig.ArgsIdxEnum.ARGS_IDX_USE_SAMPLES_IN_ORDER.idx());
        longopts[1] = new LongOpt("fail_if_updating", LongOpt.NO_ARGUMENT, null, BaseImportConfig.ArgsIdxEnum.ARGS_IDX_FAIL_IF_UPDATING.idx());
        longopts[2] = new LongOpt("interval", LongOpt.REQUIRED_ARGUMENT, null, 'L');
        longopts[3] = new LongOpt("workspace", LongOpt.REQUIRED_ARGUMENT, null, 'w');
        longopts[4] = new LongOpt("batchsize", LongOpt.REQUIRED_ARGUMENT, null, BaseImportConfig.ArgsIdxEnum.ARGS_IDX_BATCHSIZE.idx());
        longopts[5] = new LongOpt("vidmap-output", LongOpt.REQUIRED_ARGUMENT, null, BaseImportConfig.ArgsIdxEnum.ARGS_IDX_VIDMAP_OUTPUT.idx());
        longopts[6] = new LongOpt("callset-output", LongOpt.REQUIRED_ARGUMENT, null, BaseImportConfig.ArgsIdxEnum.ARGS_IDX_CALLSET_OUTPUT.idx());
        longopts[7] = new LongOpt("pass-as-bcf", LongOpt.NO_ARGUMENT, null, BaseImportConfig.ArgsIdxEnum.ARGS_IDX_PASS_AS_BCF.idx());
        longopts[8] = new LongOpt("vcf-header-output", LongOpt.REQUIRED_ARGUMENT, null, BaseImportConfig.ArgsIdxEnum.ARGS_IDX_VCF_HEADER_OUTPUT.idx());
        return longopts;
    }

    private void resolveCommandArgs(Getopt commandArgs) {
        int c;
        final int firstEnumIdx = ArgsIdxEnum.ARGS_IDX_USE_SAMPLES_IN_ORDER.idx();
        final ArgsIdxEnum[] enumArray = ArgsIdxEnum.values();
        while ((c = commandArgs.getopt()) != -1) {
            switch (c) {
                case 'w':
                    this.setWorkspace(commandArgs.getOptarg());
                    break;
                case 'L':
                    Function<String[], ChromosomeInterval> chromInterConverter = par -> new ChromosomeInterval(par[0],
                            Long.parseLong(par[1].split("-")[0]),
                            Long.parseLong(par[1].split("-")[1]));
                    this.getChromosomeIntervalList().add(chromInterConverter.apply(commandArgs.getOptarg().split(":")));
                    break;
                default: {
                    if (c >= firstEnumIdx && c < ArgsIdxEnum.ARGS_IDX_AFTER_LAST_ARG_IDX.idx()) {
                        int offset = c - firstEnumIdx;
                        assert offset < enumArray.length;
                        switch (enumArray[offset]) {
                            case ARGS_IDX_USE_SAMPLES_IN_ORDER:
                                setUseSamplesInOrderProvided(true);
                                break;
                            case ARGS_IDX_FAIL_IF_UPDATING:
                                setFailIfUpdating(true);
                                break;
                            case ARGS_IDX_BATCHSIZE:
                                setBatchSize(Integer.parseInt(commandArgs.getOptarg()));
                                break;
                            case ARGS_IDX_VIDMAP_OUTPUT:
                                setVidmapOutputFilepath(commandArgs.getOptarg());
                                break;
                            case ARGS_IDX_CALLSET_OUTPUT:
                                setCallsetOutputFilepath(commandArgs.getOptarg());
                                break;
                            case ARGS_IDX_VCF_HEADER_OUTPUT:
                                setVcfHeaderOutputFilepath(commandArgs.getOptarg());
                                break;
                            case ARGS_IDX_PASS_AS_BCF:
                                setPassAsVcf(false);
                                break;
                            default:
                                throw new IllegalArgumentException("Unknown command line option " +
                                        commandArgs.getOptarg() + " - ignored");
                        }
                    } else {
                        throw new IllegalArgumentException("Unknown command line option " +
                                commandArgs.getOptarg() + " - ignored");
                    }
                }
            }
        }
    }

    private boolean isWithinChromosomeInterval(ChromosomeInterval current, ChromosomeInterval chromInterval) {
        return (current.getStart() >= chromInterval.getStart() &&
                current.getStart() <= chromInterval.getEnd() &&
                current.getContig().equals(chromInterval.getContig())) ||
                (current.getEnd() >= chromInterval.getStart() &&
                        current.getEnd() <= chromInterval.getEnd() &&
                        current.getContig().equals(chromInterval.getContig()));
    }

    private boolean isThereChromosomeIntervalIntersection(List<ChromosomeInterval> chromIntervals, boolean isThereInter) {
        if (chromIntervals.isEmpty() || chromIntervals.size() < 2) return isThereInter;
        ChromosomeInterval head = chromIntervals.get(0);
        List<ChromosomeInterval> tail = chromIntervals.subList(1, chromIntervals.size());

        for (ChromosomeInterval chrom : tail) {
            boolean interEval = isWithinChromosomeInterval(head, chrom);
            isThereInter = isThereInter || interEval;
        }

        return isThereChromosomeIntervalIntersection(tail, isThereInter);
    }

    private boolean isThereChromosomeIntervalIntersection(List<ChromosomeInterval> chromIntervals) {
        return isThereChromosomeIntervalIntersection(chromIntervals, false);
    }

    public List<ChromosomeInterval> getChromosomeIntervalList() {
        return chromosomeIntervalList;
    }

    public String getWorkspace() {
        return workspace;
    }

    public Long getVcfBufferSizePerColumnPartition() {
        return vcfBufferSizePerColumnPartition;
    }

    public Long getSegmentSize() {
        return segmentSize;
    }

    public Long getLbRowIdx() { return lbRowIdx; }

    public Long getUbRowIdx() {
        return ubRowIdx;
    }

    public boolean isUseSamplesInOrderProvided() {
        return useSamplesInOrderProvided;
    }

    public boolean isFailIfUpdating() {
        return failIfUpdating;
    }

    public int getRank() {
        return rank;
    }

    public boolean isValidateSampleToReaderMap() {
        return validateSampleToReaderMap;
    }

    public boolean isPassAsVcf() {
        return passAsVcf;
    }

    public int getBatchSize() {
        return batchSize;
    }

    public String getVidmapOutputFilepath() {
        return vidmapOutputFilepath;
    }

    public String getCallsetOutputFilepath() {
        return callsetOutputFilepath;
    }

    public String getVcfHeaderOutputFilepath() {
        return vcfHeaderOutputFilepath;
    }

    public List<String> getFiles() { return files; }

    public void setChromosomeIntervalList(List<ChromosomeInterval> chromosomeIntervalList) {
        this.chromosomeIntervalList = chromosomeIntervalList;
    }

    public void setWorkspace(String workspace) {
        this.workspace = workspace;
    }

    public void setVcfBufferSizePerColumnPartition(Long vcfBufferSizePerColumnPartition) {
        this.vcfBufferSizePerColumnPartition = vcfBufferSizePerColumnPartition;
    }

    public void setSegmentSize(Long segmentSize) {
        this.segmentSize = segmentSize;
    }

    public void setLbRowIdx(Long lbRowIdx) {
        this.lbRowIdx = lbRowIdx;
    }

    public void setUbRowIdx(Long ubRowIdx) {
        this.ubRowIdx = ubRowIdx;
    }

    public void setUseSamplesInOrderProvided(boolean useSamplesInOrderProvided) {
        this.useSamplesInOrderProvided = useSamplesInOrderProvided;
    }

    public void setFailIfUpdating(boolean failIfUpdating) {
        this.failIfUpdating = failIfUpdating;
    }

    public void setRank(int rank) {
        this.rank = rank;
    }

    public void setValidateSampleToReaderMap(boolean validateSampleToReaderMap) {
        this.validateSampleToReaderMap = validateSampleToReaderMap;
    }

    public void setPassAsVcf(boolean passAsVcf) {
        this.passAsVcf = passAsVcf;
    }

    public void setBatchSize(int batchSize) {
        this.batchSize = batchSize;
    }

    public void setVidmapOutputFilepath(String vidmapOutputFilepath) {
        this.vidmapOutputFilepath = vidmapOutputFilepath;
    }

    public void setCallsetOutputFilepath(String callsetOutputFilepath) {
        this.callsetOutputFilepath = callsetOutputFilepath;
    }

    public void setVcfHeaderOutputFilepath(String vcfHeaderOutputFilepath) {
        this.vcfHeaderOutputFilepath = vcfHeaderOutputFilepath;
    }

    public void setFiles(List<String> files) {
        this.files = files;
    }

    public enum ArgsIdxEnum {
        ARGS_IDX_USE_SAMPLES_IN_ORDER(1000),
        ARGS_IDX_FAIL_IF_UPDATING(1001),
        ARGS_IDX_BATCHSIZE(1002),
        ARGS_IDX_VIDMAP_OUTPUT(1003),
        ARGS_IDX_CALLSET_OUTPUT(1004),
        ARGS_IDX_PASS_AS_BCF(1005),
        ARGS_IDX_VCF_HEADER_OUTPUT(1006),
        ARGS_IDX_AFTER_LAST_ARG_IDX(1007);

        private final int mArgsIdx;

        ArgsIdxEnum(final int idx) {
            mArgsIdx = idx;
        }

        public int idx() {
            return mArgsIdx;
        }
    }

}
