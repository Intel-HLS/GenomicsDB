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

import com.intel.genomicsdb.importer.model.ChromosomeInterval;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.nio.file.Path;
import java.util.*;

public class ParallelImportConfig {
    private GenomicsDBImportConfiguration.ImportConfiguration importConfiguration;
    private List<ChromosomeInterval> chromosomeIntervalList = new ArrayList<>();
    private int rank = 0;
    private boolean validateSampleToReaderMap;
    private boolean passAsVcf = true;
    private int batchSize = 1000000;
    private Set<VCFHeaderLine> mergedHeader;
    private Map<String, Path> sampleNameToVcfPath;
    private Func<Map<String, Path>, Integer, Integer, Map<String, FeatureReader<VariantContext>>> sampleToReaderMap;

    public ParallelImportConfig(final GenomicsDBImportConfiguration.ImportConfiguration importConfiguration,
                                final List<ChromosomeInterval> chromosomeIntervalList,
                                final int rank,
                                final boolean validateSampleToReaderMap,
                                final boolean passAsVcf,
                                final int batchSize,
                                final Set<VCFHeaderLine> mergedHeader,
                                final Map<String, Path> sampleNameToVcfPath,
                                final Func<Map<String, Path>, Integer, Integer,
                                    Map<String, FeatureReader<VariantContext>>> sampleToReaderMap) {
        this.setImportConfiguration(importConfiguration);
        this.setChromosomeIntervalList(chromosomeIntervalList);
        this.setRank(rank);
        this.setValidateSampleToReaderMap(validateSampleToReaderMap);
        this.setPassAsVcf(passAsVcf);
        this.setBatchSize(batchSize);
        this.setMergedHeader(mergedHeader);
        this.setSampleNameToVcfPath(sampleNameToVcfPath);
        this.setSampleToReaderMap(sampleToReaderMap);
    }

    protected ParallelImportConfig() {
    }

    private boolean isWithinChromosomeInterval(final ChromosomeInterval current, final ChromosomeInterval chromInterval) {
        return (current.getStart() >= chromInterval.getStart() &&
                current.getStart() <= chromInterval.getEnd() &&
                current.getContig().equals(chromInterval.getContig())) ||
                (current.getEnd() >= chromInterval.getStart() &&
                        current.getEnd() <= chromInterval.getEnd() &&
                        current.getContig().equals(chromInterval.getContig()));
    }

    private boolean isThereChromosomeIntervalIntersection(final List<ChromosomeInterval> chromIntervals, boolean isThereInter) {
        if (chromIntervals.isEmpty() || chromIntervals.size() < 2) return isThereInter;
        ChromosomeInterval head = chromIntervals.get(0);
        List<ChromosomeInterval> tail = chromIntervals.subList(1, chromIntervals.size());

        for (ChromosomeInterval chrom : tail) {
            boolean interEval = isWithinChromosomeInterval(head, chrom);
            isThereInter = isThereInter || interEval;
        }

        return isThereChromosomeIntervalIntersection(tail, isThereInter);
    }

    private boolean isThereChromosomeIntervalIntersection(final List<ChromosomeInterval> chromIntervals) {
        return isThereChromosomeIntervalIntersection(chromIntervals, false);
    }

    void validateChromosomeIntervals(final List<ChromosomeInterval> chromosomeIntervalList) {
        if (isThereChromosomeIntervalIntersection(chromosomeIntervalList))
            throw new IllegalArgumentException("There are multiple intervals sharing same value. This is not allowed. " +
                    "Intervals should be defined without intersections.");
    }

    public GenomicsDBImportConfiguration.ImportConfiguration getImportConfiguration() {
        return importConfiguration;
    }

    public void setImportConfiguration(GenomicsDBImportConfiguration.ImportConfiguration importConfiguration) {
        this.importConfiguration = importConfiguration;
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

    public Func<Map<String, Path>, Integer, Integer, Map<String, FeatureReader<VariantContext>>> sampleToReaderMapCreator() {
        return sampleToReaderMap;
    }

    public void setSampleToReaderMap(Func<Map<String, Path>, Integer, Integer, Map<String, FeatureReader<VariantContext>>> sampleToReaderMap) {
        this.sampleToReaderMap = sampleToReaderMap;
    }

    public List<ChromosomeInterval> getChromosomeIntervalList() {
        return chromosomeIntervalList;
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

    public void setChromosomeIntervalList(List<ChromosomeInterval> chromosomeIntervalList) {
        this.validateChromosomeIntervals(chromosomeIntervalList);
        this.chromosomeIntervalList = chromosomeIntervalList;
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

    public Set<VCFHeaderLine> getMergedHeader() {
        return mergedHeader;
    }

    public void setMergedHeader(Set<VCFHeaderLine> mergedHeader) {
        this.mergedHeader = mergedHeader;
    }

}
