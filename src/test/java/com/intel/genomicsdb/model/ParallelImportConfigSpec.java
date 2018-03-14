package com.intel.genomicsdb.model;

import com.intel.genomicsdb.importer.model.ChromosomeInterval;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.util.*;

public class ParallelImportConfigSpec {

    @Test(testName = "should throw an exception when there is an intersection between chromosome intervals",
            expectedExceptions = IllegalArgumentException.class,
            expectedExceptionsMessageRegExp = "There are multiple intervals sharing same value. This is not allowed. " +
                    "Intervals should be defined without intersections.")
    public void shouldThrowExceptionWhenThereIsAnIntersectionBetweenChromosomeIntervals() {
        //Given
        ChromosomeInterval chromosomeInterval1 = new ChromosomeInterval("1", 1, 10);
        ChromosomeInterval chromosomeInterval2 = new ChromosomeInterval("1", 5, 15);
        final List<ChromosomeInterval> chromosomeIntervals = new ArrayList<>();
        chromosomeIntervals.add(chromosomeInterval1);
        chromosomeIntervals.add(chromosomeInterval2);

        //When
        createBaseImportConfig(chromosomeIntervals);

        //Then
        //Exception is expected
    }

    @Test(testName = "should throw an exception when there are not intersections between chromosome intervals")
    public void shouldNotThrowExceptionWhenThereAreIntersectionsButDifferentChromosomes() {
        //Given
        ChromosomeInterval chromosomeInterval1 = new ChromosomeInterval("1", 1, 10);
        ChromosomeInterval chromosomeInterval2 = new ChromosomeInterval("2", 5, 15);
        final List<ChromosomeInterval> chromosomeIntervals = new ArrayList<>();
        chromosomeIntervals.add(chromosomeInterval1);
        chromosomeIntervals.add(chromosomeInterval2);

        //When
        ParallelImportConfig config = createBaseImportConfig(chromosomeIntervals);

        //Then
        Assert.assertEquals(config.getChromosomeIntervalList().size(), 2);
    }

    private ParallelImportConfig createBaseImportConfig(final List<ChromosomeInterval> chromosomeIntervalList) {
        List<ChromosomeInterval> chromosomeIntervals = chromosomeIntervalList;
        GenomicsDBImportConfiguration.ImportConfiguration configuration = GenomicsDBImportConfiguration.ImportConfiguration.getDefaultInstance();
        int rank = 1;
        boolean validateSampleToReaderMap = true;
        boolean passAsVcf = true;
        int batchSize = 1000;
        Set<VCFHeaderLine> mergedHeader = new HashSet();
        Map<String, Path> sampleNameToVcfPath = new TreeMap<>();
        ParallelImportConfig.Func<Map<String, Path>, Integer, Integer,
                Map<String, FeatureReader<VariantContext>>> sampleToReaderMap = (a, b, c) -> new TreeMap<>();

        return new ParallelImportConfig(configuration, chromosomeIntervals, rank, validateSampleToReaderMap,
                passAsVcf, batchSize, mergedHeader, sampleNameToVcfPath, sampleToReaderMap);
    }
}
