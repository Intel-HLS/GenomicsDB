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
import java.util.function.Function;
import java.lang.Integer;
import java.util.stream.IntStream;
import java.util.stream.Collectors;

/**
 * This implementation extends what is in GenomicsDBImportConfiguration. Add extra data that is needed for parallel
 * import.
 */
public class ImportConfig {
    private GenomicsDBImportConfiguration.ImportConfiguration importConfiguration;
    private boolean validateSampleToReaderMap = false;
    private boolean passAsVcf = true;
    private boolean useSamplesInOrder = false;
    private int batchSize = 0;
    private Set<VCFHeaderLine> mergedHeader;
    private Map<String, Path> sampleNameToVcfPath;
    private Func<Map<String, Path>, Integer, Integer, Map<String, FeatureReader<VariantContext>>> sampleToReaderMapCreator;
    private Function<BatchCompletionCallbackFunctionArgument, Void> functionToCallOnBatchCompletion = null;
    private String outputVidMapJsonFile = null;
    private String outputCallsetMapJsonFile = null;
    private String outputVcfHeaderFile = null;

    /**
     * Main ImportConfig constructor
     *
     * @param importConfiguration       GenomicsDBImportConfiguration protobuf object
     * @param validateSampleToReaderMap Flag for validating sample to reader map
     * @param passAsVcf                 Flag for indicating that a VCF is being passed
     * @param batchSize                 Batch size
     * @param mergedHeader              Required header
     * @param sampleNameToVcfPath       Sample name to VCF path map
     * @param sampleToReaderMapCreator  Function used for creating sampleToReaderMap
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

    class SortGenomicsDBPartition implements Comparator<Integer> {

        public SortGenomicsDBPartition(final List<GenomicsDBImportConfiguration.Partition> partitionList) {
            this.partitionList = partitionList;
        }

        public int compare(final Integer lIdx, final Integer rIdx) {
            final GenomicsDBImportConfiguration.Partition l = partitionList.get(lIdx);
            final GenomicsDBImportConfiguration.Partition r = partitionList.get(rIdx);
            int contigCompare = l.getBegin().getContigPosition().getContig().compareTo(r.getBegin().getContigPosition().getContig());
            if(contigCompare < 0)
                return -1;
            if(contigCompare > 0)
                return 1;
            long lPos = l.getBegin().getContigPosition().getPosition();
            long rPos = r.getBegin().getContigPosition().getPosition();
            return (lPos < rPos) ? -1 : (lPos > rPos) ? 1 : 0;
        }

        private List<GenomicsDBImportConfiguration.Partition> partitionList = null;
    }

    private boolean isThereChromosomeIntervalIntersection(final List<Integer> sortedPartitionIdxs) {
        List<GenomicsDBImportConfiguration.Partition> partitions = this.importConfiguration.getColumnPartitionsList();
        for(int i=0;i<sortedPartitionIdxs.size()-1;++i) {
            Coordinates.ContigPosition currBegin = partitions.get(sortedPartitionIdxs.get(i)).getBegin().getContigPosition();
            Coordinates.ContigPosition currEnd = partitions.get(sortedPartitionIdxs.get(i)).hasEnd()
                ? partitions.get(sortedPartitionIdxs.get(i)).getEnd().getContigPosition()
                : null;
            Coordinates.ContigPosition nextBegin = partitions.get(sortedPartitionIdxs.get(i+1)).getBegin().getContigPosition();
            if(currBegin.getContig().equals(nextBegin.getContig()) //same contig
                    && (currEnd == null ||      //the first interval spans the full contig
                        currEnd.getPosition() >= nextBegin.getPosition())) //overlap
                return true;
        }
        return false;
    }

    private boolean isThereChromosomeIntervalIntersection() {
        List<GenomicsDBImportConfiguration.Partition> partitions = this.importConfiguration.getColumnPartitionsList();
        List<Integer> partitionIdxList = IntStream.range(0, partitions.size()).boxed().collect(Collectors.toList());
        Collections.sort(partitionIdxList, new SortGenomicsDBPartition(partitions));
        return isThereChromosomeIntervalIntersection(partitionIdxList);
    }

    void validateChromosomeIntervals() {
        for(GenomicsDBImportConfiguration.Partition currPartition : importConfiguration.getColumnPartitionsList()) {
            if(!currPartition.getBegin().hasContigPosition() || (currPartition.hasEnd() && !currPartition.getEnd().hasContigPosition()))
                throw new IllegalArgumentException("Must use contig positions while using multi-interval import");
            if(currPartition.hasEnd()) {
               if(!currPartition.getBegin().getContigPosition().getContig().equals(
                           currPartition.getEnd().getContigPosition().getContig()))
                   throw new IllegalArgumentException("Both begin and end for a partition must be in the same contig");
               if(currPartition.getBegin().getContigPosition().getPosition() >
                       currPartition.getEnd().getContigPosition().getPosition())
                   throw new IllegalArgumentException("End of a partition cannot be less than begin");
            }
        }
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

    public String getOutputVcfHeaderFile() {
        return outputVcfHeaderFile;
    }

    public void setOutputVcfHeaderFile(String outputVcfHeaderFile) {
        this.outputVcfHeaderFile = outputVcfHeaderFile;
    }

    @FunctionalInterface
    public interface Func<T1, T2, T3, R> {
        R apply(T1 t1, T2 t2, T3 t3);
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

    public void setFunctionToCallOnBatchCompletion(Function<BatchCompletionCallbackFunctionArgument, Void> functionToCallOnBatchCompletion) {
        this.functionToCallOnBatchCompletion = functionToCallOnBatchCompletion;
    }

    public Function<BatchCompletionCallbackFunctionArgument, Void> getFunctionToCallOnBatchCompletion() {
        return functionToCallOnBatchCompletion;
    }
}
