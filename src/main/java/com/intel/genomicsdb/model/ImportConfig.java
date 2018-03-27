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

import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.nio.file.Path;
import java.util.*;
import java.lang.Integer;

/**
 * This implementation extends what is in GenomicsDBImportConfiguration. Add extra data that is needed for parallel
 * import.
 */
public class ImportConfig {
    private GenomicsDBImportConfiguration.ImportConfiguration importConfiguration;
    private boolean validateSampleToReaderMap = false;
    private boolean passAsVcf = true;
    private boolean useSamplesInOrder = false;
    private int batchSize = Integer.MAX_VALUE;
    private Set<VCFHeaderLine> mergedHeader;
    private Map<String, Path> sampleNameToVcfPath;
    private Func<Map<String, Path>, Integer, Integer, Map<String, FeatureReader<VariantContext>>> sampleToReaderMapCreator;
    private String outputVidMapJsonFile = null;
    private String outputCallsetMapJsonFile = null;

    /**
     * Main ImportConfig constructor
     * @param importConfiguration GenomicsDBImportConfiguration protobuf object
     * @param validateSampleToReaderMap Flag for validating sample to reader map
     * @param passAsVcf Flag for indicating that a VCF is being passed
     * @param batchSize Batch size
     * @param mergedHeader Required header
     * @param sampleNameToVcfPath Sample name to VCF path map
     * @param sampleToReaderMapCreator Function used for creating sampleToReaderMap
     */
    public ImportConfig(final GenomicsDBImportConfiguration.ImportConfiguration importConfiguration,
                        final boolean validateSampleToReaderMap,
                        final boolean passAsVcf,
                        final int batchSize,
                        final Set<VCFHeaderLine> mergedHeader,
                        final Map<String, Path> sampleNameToVcfPath,
                        final Func<Map<String, Path>, Integer, Integer,
                                    Map<String, FeatureReader<VariantContext>>> sampleToReaderMapCreator) {
        this.setImportConfiguration(importConfiguration);
        this.validateChromosomeIntervals();
        this.setValidateSampleToReaderMap(validateSampleToReaderMap);
        this.setPassAsVcf(passAsVcf);
        this.setBatchSize(batchSize);
        this.setMergedHeader(mergedHeader);
        this.setSampleNameToVcfPath(sampleNameToVcfPath);
        this.setSampleToReaderMapCreator(sampleToReaderMapCreator);
    }

    /**
     * Deep copy constructor
     * @param source Source ImportConfig
     */
    public ImportConfig(final ImportConfig source) {
        this.setImportConfiguration(source.getImportConfiguration() != null ?
                source.getImportConfiguration().toBuilder().build() : null);
        this.validateChromosomeIntervals();
        this.setValidateSampleToReaderMap(source.isValidateSampleToReaderMap());
        this.setPassAsVcf(source.isPassAsVcf());
        this.setUseSamplesInOrder(source.isUseSamplesInOrder());
        this.setBatchSize(source.getBatchSize());
        //Cannot deep copy mergedHeader easily - so do shallow copy and make it unmodifiable
        this.setMergedHeader(source.getMergedHeader());
        //Deep copy sample name to path
        if(source.getSampleNameToVcfPath() != null) {
            this.setSampleNameToVcfPath(new LinkedHashMap<>(source.getSampleNameToVcfPath()));
        }
        //Function object - no deep copy
        this.setSampleToReaderMapCreator(source.sampleToReaderMapCreator());
        //Deep copy
        this.setOutputVidmapJsonFile(source.getOutputVidmapJsonFile() != null ?
                new String(source.getOutputVidmapJsonFile()) : null);
        this.setOutputCallsetmapJsonFile(source.getOutputCallsetmapJsonFile() != null ?
                new String(source.getOutputCallsetmapJsonFile()) : null);
    }

    protected ImportConfig() {
    }

    private boolean isWithinChromosomeInterval(final GenomicsDBImportConfiguration.Partition current,
                                               final GenomicsDBImportConfiguration.Partition chromInterval) {
        return (current.getBegin().getContigPosition().getPosition() >= chromInterval.getBegin().getContigPosition().getPosition() &&
                current.getBegin().getContigPosition().getPosition() <= chromInterval.getEnd().getContigPosition().getPosition() &&
                current.getBegin().getContigPosition().getContig().equals(chromInterval.getBegin().getContigPosition().getContig())) ||
                (current.getEnd().getContigPosition().getPosition() >= chromInterval.getBegin().getContigPosition().getPosition() &&
                        current.getEnd().getContigPosition().getPosition() <= chromInterval.getEnd().getContigPosition().getPosition() &&
                        current.getBegin().getContigPosition().getContig().equals(chromInterval.getBegin().getContigPosition().getContig()));
    }

    private boolean isThereChromosomeIntervalIntersection(final List<GenomicsDBImportConfiguration.Partition> chromIntervals, boolean isThereInter) {
        if (chromIntervals.isEmpty() || chromIntervals.size() < 2) return isThereInter;
        GenomicsDBImportConfiguration.Partition head = chromIntervals.get(0);
        List<GenomicsDBImportConfiguration.Partition> tail = chromIntervals.subList(1, chromIntervals.size());

        for (GenomicsDBImportConfiguration.Partition chrom : tail) {
            boolean interEval = isWithinChromosomeInterval(head, chrom);
            isThereInter = isThereInter || interEval;
        }

        return isThereChromosomeIntervalIntersection(tail, isThereInter);
    }

    private boolean isThereChromosomeIntervalIntersection() {
        List<GenomicsDBImportConfiguration.Partition> partitions = this.importConfiguration.getColumnPartitionsList();
        return isThereChromosomeIntervalIntersection(partitions, false);
    }

    void validateChromosomeIntervals() {
        if (isThereChromosomeIntervalIntersection())
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
        return sampleToReaderMapCreator;
    }

    public void setSampleToReaderMapCreator(
            Func<Map<String, Path>, Integer, Integer, Map<String, FeatureReader<VariantContext>>> sampleToReaderMapCreator) {
        this.sampleToReaderMapCreator = sampleToReaderMapCreator;
    }

    public boolean isValidateSampleToReaderMap() {
        return validateSampleToReaderMap;
    }

    public boolean isPassAsVcf() {
        return passAsVcf;
    }

    public boolean isUseSamplesInOrder() {
        return useSamplesInOrder;
    }

    public int getBatchSize() {
        return batchSize;
    }

    public void setValidateSampleToReaderMap(boolean validateSampleToReaderMap) {
        this.validateSampleToReaderMap = validateSampleToReaderMap;
    }

    public void setPassAsVcf(boolean passAsVcf) {
        this.passAsVcf = passAsVcf;
    }

    public void setUseSamplesInOrder(final boolean useSamplesInOrder) {
        this.useSamplesInOrder = useSamplesInOrder;
    }

    public void setBatchSize(int batchSize) {
        this.batchSize = batchSize;
    }

    public Set<VCFHeaderLine> getMergedHeader() {
        return Collections.unmodifiableSet(mergedHeader);
    }

    public void setMergedHeader(Set<VCFHeaderLine> mergedHeader) {
        this.mergedHeader = mergedHeader;
    }

    public String getOutputCallsetmapJsonFile() {
        return outputCallsetMapJsonFile;
    }

    public String getOutputVidmapJsonFile() {
        return outputVidMapJsonFile;
    }

    public void setOutputCallsetmapJsonFile(final String outputCallsetMapJsonFile) {
        this.outputCallsetMapJsonFile = outputCallsetMapJsonFile;
    }

    public void setOutputVidmapJsonFile(final String outputVidMapJsonFile) {
        this.outputVidMapJsonFile = outputVidMapJsonFile;
    }
}
