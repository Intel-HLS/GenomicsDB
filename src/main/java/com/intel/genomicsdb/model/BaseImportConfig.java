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
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;

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
    private Set<VCFHeaderLine> mergedHeader;
    private Map<String, Path> sampleNameToVcfPath;
    private Func<Map<String, Path>, Integer, Integer, Map<String, FeatureReader<VariantContext>>> sampleToReaderMap;

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
                            final Set<VCFHeaderLine> mergedHeader,
                            final Map<String, Path> sampleNameToVcfPath,
                            final Func<Map<String, Path>, Integer, Integer,
                                    Map<String, FeatureReader<VariantContext>>> sampleToReaderMap) {
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
        this.setMergedHeader(mergedHeader);
        this.setSampleNameToVcfPath(sampleNameToVcfPath);
        this.setSampleToReaderMap(sampleToReaderMap);
    }

    protected BaseImportConfig() {
    }

    private <T> void validate(T attribute, Function<T, Boolean> validation, String errorMessage) {
        if(!validation.apply(attribute)) throw new IllegalArgumentException(errorMessage);
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

    void validateChromosomeIntervals(List<ChromosomeInterval> chromosomeIntervalList) {
        if (isThereChromosomeIntervalIntersection(chromosomeIntervalList))
            throw new IllegalArgumentException("There are multiple intervals sharing same value. This is not allowed. " +
                    "Intervals should be defined without intersections.");
    }

    public boolean isUseSamplesInOrderProvided() {
        return useSamplesInOrderProvided;
    }

    public void setUseSamplesInOrderProvided(boolean useSamplesInOrderProvided) {
        this.useSamplesInOrderProvided = useSamplesInOrderProvided;
    }

    @FunctionalInterface
    public interface Func<T1, T2, T3, R>{
        R apply(T1 t1,T2 t2,T3 t3);
    }

    public Map<String, Path> getSampleNameToVcfPath() {
        return sampleNameToVcfPath;
    }

    public void setSampleNameToVcfPath(Map<String, Path> sampleNameToVcfPath) {
        this.sampleNameToVcfPath = sampleNameToVcfPath;
    }

    public Func<Map<String, Path>, Integer, Integer, Map<String, FeatureReader<VariantContext>>> getSampleToReaderMap() {
        return sampleToReaderMap;
    }

    public void setSampleToReaderMap(Func<Map<String, Path>, Integer, Integer, Map<String, FeatureReader<VariantContext>>> sampleToReaderMap) {
        this.sampleToReaderMap = sampleToReaderMap;
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

    public void setChromosomeIntervalList(List<ChromosomeInterval> chromosomeIntervalList) {
        this.validateChromosomeIntervals(chromosomeIntervalList);
        this.chromosomeIntervalList = chromosomeIntervalList;
    }

    public void setWorkspace(String workspace) {
        this.validate(workspace, (value) -> value != null && !value.isEmpty(), "Workspace name cannot be null or empty.");
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

    public Set<VCFHeaderLine> getMergedHeader() {
        return mergedHeader;
    }

    public void setMergedHeader(Set<VCFHeaderLine> mergedHeader) {
        this.mergedHeader = mergedHeader;
    }

}
