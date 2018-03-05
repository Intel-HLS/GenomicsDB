package com.intel.genomicsdb.model;

import com.intel.genomicsdb.ChromosomeInterval;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.util.*;

public class BaseImportConfigSpec {

    @Test(testName = "should throw an exception when workspace is null or empty",
            expectedExceptions = IllegalArgumentException.class,
            expectedExceptionsMessageRegExp = "Workspace name cannot be null or empty.")
    public void shouldThrowExceptionWhenWorkspaceIsNullOrEmpty() {
        //Given
        String workspace = "";
        ChromosomeInterval chromosomeInterval1 = new ChromosomeInterval("1", 1, 10);
        ChromosomeInterval chromosomeInterval2 = new ChromosomeInterval("1", 11, 20);
        final List<ChromosomeInterval> chromosomeIntervals = new ArrayList<>();
        chromosomeIntervals.add(chromosomeInterval1);
        chromosomeIntervals.add(chromosomeInterval2);

        //When
        createBaseImportConfig(chromosomeIntervals, workspace);

        //Then
        //Exception is expected
    }

    @Test(testName = "should throw an exception when there is an intersection between chromosome intervals",
            expectedExceptions = IllegalArgumentException.class,
            expectedExceptionsMessageRegExp = "There are multiple intervals sharing same value. This is not allowed. " +
                    "Intervals should be defined without intersections.")
    public void shouldThrowExceptionWhenThereIsAnIntersectionBetweenChromosomeIntervals() {
        //Given
        String workspace = "path/to/ws";
        ChromosomeInterval chromosomeInterval1 = new ChromosomeInterval("1", 1, 10);
        ChromosomeInterval chromosomeInterval2 = new ChromosomeInterval("1", 5, 15);
        final List<ChromosomeInterval> chromosomeIntervals = new ArrayList<>();
        chromosomeIntervals.add(chromosomeInterval1);
        chromosomeIntervals.add(chromosomeInterval2);

        //When
        createBaseImportConfig(chromosomeIntervals, workspace);

        //Then
        //Exception is expected
    }

    @Test(testName = "should throw an exception when there are not intersections between chromosome intervals")
    public void shouldNotThrowExceptionWhenThereAreIntersectionsButDifferentChromosomes() {
        //Given
        String workspace = "path/to/ws";
        ChromosomeInterval chromosomeInterval1 = new ChromosomeInterval("1", 1, 10);
        ChromosomeInterval chromosomeInterval2 = new ChromosomeInterval("2", 5, 15);
        final List<ChromosomeInterval> chromosomeIntervals = new ArrayList<>();
        chromosomeIntervals.add(chromosomeInterval1);
        chromosomeIntervals.add(chromosomeInterval2);

        //When
        BaseImportConfig config = createBaseImportConfig(chromosomeIntervals, workspace);

        //Then
        Assert.assertEquals(config.getChromosomeIntervalList().size(), 2);
    }

    private BaseImportConfig createBaseImportConfig(final List<ChromosomeInterval> chromosomeIntervalList,
                                                    final String ws) {
        String workspace = ws;
        List<ChromosomeInterval> chromosomeIntervals = chromosomeIntervalList;
        long vcfBufferSizePerColumnPartition = 1;
        long segmentSize = 1;
        long lbRowIdx = 1;
        long ubRowIdx = 2;
        boolean useSamplesInOrderProvided = false;
        boolean failIfUpdating = false;
        int rank = 1;
        boolean validateSampleToReaderMap = true;
        boolean passAsVcf = true;
        int batchSize = 1000;
        String vidmapOutputFilepath = "/path/to/somewhere";
        String callsetOutputFilepath = "/path/to/somewhere";
        String vcfHeaderOutputFilepath = "/path/to/somewhere";
        Set<VCFHeaderLine> mergedHeader = new HashSet();
        Map<String, Path> sampleNameToVcfPath = new TreeMap<>();
        BaseImportConfig.Func<Map<String, Path>, Integer, Integer,
                Map<String, FeatureReader<VariantContext>>> sampleToReaderMap = (a, b, c) -> new TreeMap<>();

        return new BaseImportConfig(chromosomeIntervals, workspace, vcfBufferSizePerColumnPartition, segmentSize,
                lbRowIdx, ubRowIdx, useSamplesInOrderProvided, failIfUpdating, rank, validateSampleToReaderMap,
                passAsVcf, batchSize, vidmapOutputFilepath, callsetOutputFilepath, vcfHeaderOutputFilepath,
                mergedHeader, sampleNameToVcfPath, sampleToReaderMap);
    }
}
