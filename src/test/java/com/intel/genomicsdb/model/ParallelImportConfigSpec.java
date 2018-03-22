package com.intel.genomicsdb.model;

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
        //When
        createBaseImportConfig(true);

        //Then
        //Exception is expected
    }

    @Test(testName = "should not throw an exception when there are not intersections between chromosome intervals")
    public void shouldNotThrowExceptionWhenThereAreIntersectionsButDifferentChromosomes() {
        //Given
        //When
        ParallelImportConfig config = createBaseImportConfig(false);

        //Then
        Assert.assertEquals(config.getImportConfiguration().getColumnPartitionsList().size(), 2);
    }

    private ParallelImportConfig createBaseImportConfig(final boolean withIntersection) {
        String defaultName = "a";
        String secondPartition = withIntersection ? defaultName : "b";
        GenomicsDBImportConfiguration.ImportConfiguration configuration = GenomicsDBImportConfiguration.ImportConfiguration.newBuilder()
                .addColumnPartitions(GenomicsDBImportConfiguration.Partition.newBuilder().setBegin(
                        Coordinates.GenomicsDBColumn.newBuilder().setContigPosition(
                                Coordinates.ContigPosition.newBuilder().setContig("a").setPosition(1).build()
                        ).build()
                        ).setEnd(
                        Coordinates.GenomicsDBColumn.newBuilder().setContigPosition(
                                Coordinates.ContigPosition.newBuilder().setContig("a").setPosition(2).build()
                        ).build()
                        ).build()
                )
                .addColumnPartitions(GenomicsDBImportConfiguration.Partition.newBuilder().setBegin(
                        Coordinates.GenomicsDBColumn.newBuilder().setContigPosition(
                                Coordinates.ContigPosition.newBuilder().setContig(secondPartition).setPosition(1).build()
                        ).build()
                        ).setEnd(
                        Coordinates.GenomicsDBColumn.newBuilder().setContigPosition(
                                Coordinates.ContigPosition.newBuilder().setContig(secondPartition).setPosition(2).build()
                        ).build()
                        ).build()
                )
                .setSizePerColumnPartition(16000).build();
        boolean validateSampleToReaderMap = true;
        boolean passAsVcf = true;
        int batchSize = 1000;
        Set<VCFHeaderLine> mergedHeader = new HashSet();
        Map<String, Path> sampleNameToVcfPath = new TreeMap<>();
        ParallelImportConfig.Func<Map<String, Path>, Integer, Integer,
                Map<String, FeatureReader<VariantContext>>> sampleToReaderMap = (a, b, c) -> new TreeMap<>();

        return new ParallelImportConfig(configuration, validateSampleToReaderMap, passAsVcf, batchSize, mergedHeader,
                sampleNameToVcfPath, sampleToReaderMap);
    }
}
