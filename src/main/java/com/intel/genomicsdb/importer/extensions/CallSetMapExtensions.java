package com.intel.genomicsdb.importer.extensions;

import com.intel.genomicsdb.model.GenomicsDBCallsetsMapProto;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public interface CallSetMapExtensions {
    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes.
     *
     * @param sampleNames               list of sample names
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @return Mappings of callset (sample) names to TileDB rows
     */
    default GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(final List<String> sampleNames,
                                                                                 final boolean useSamplesInOrderProvided) {
        return generateSortedCallSetMap(sampleNames, useSamplesInOrderProvided, 0L);
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes. This method is synchronized
     * to block multiple invocations (if by any chance) disturb the order in
     * which TileDB row indexes are generated
     *
     * @param sampleNames               list of sample names
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @param lbRowIdx                  Smallest row idx which should be imported by this object
     * @return Mappings of callset (sample) names to TileDB rows
     */
    default GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(final List<String> sampleNames,
                                                                                 final boolean useSamplesInOrderProvided,
                                                                                 final long lbRowIdx) {
        if (!useSamplesInOrderProvided) Collections.sort(sampleNames);

        GenomicsDBCallsetsMapProto.CallsetMappingPB.Builder callsetMapBuilder =
                GenomicsDBCallsetsMapProto.CallsetMappingPB.newBuilder();

        long tileDBRowIndex = lbRowIdx;

        for (String sampleName : sampleNames) {
            GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap.Builder idMapBuilder =
                    GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap.newBuilder();

            idMapBuilder.setSampleName(sampleName).setRowIdx(tileDBRowIndex++).setIdxInFile(0)
                    .setStreamName(sampleName + "_stream");

            GenomicsDBCallsetsMapProto.SampleIDToTileDBIDMap sampleIDToTileDBIDMap =
                    idMapBuilder.build();

            callsetMapBuilder.addCallsets(sampleIDToTileDBIDMap);
        }

        return callsetMapBuilder.build();
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes.
     * <p>
     * Assume one sample per input GVCF file
     *
     * @param sampleToReaderMap         Variant Readers objects of the input GVCF files
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @param validateSampleToReaderMap
     * @return Mappings of callset (sample) names to TileDB rows
     */
    default GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(
            final Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
            final boolean validateSampleToReaderMap,
            final boolean useSamplesInOrderProvided) {
        return generateSortedCallSetMap(sampleToReaderMap, useSamplesInOrderProvided, validateSampleToReaderMap, 0L);
    }

    /**
     * Creates a sorted list of callsets and generates unique TileDB
     * row indices for them. Sorted to maintain order between
     * distributed share-nothing load processes.
     * <p>
     * Assume one sample per input GVCF file
     *
     * @param sampleToReaderMap         Variant Readers objects of the input GVCF files
     * @param useSamplesInOrderProvided If True, do not sort the samples,
     *                                  use the order they appear in
     * @param validateSampleMap         Check i) whether sample names are consistent
     *                                  with headers and ii) feature readers are valid
     *                                  in sampleToReaderMap
     * @param lbRowIdx                  Smallest row idx which should be imported by this object
     * @return Mappings of callset (sample) names to TileDB rows
     */
    default GenomicsDBCallsetsMapProto.CallsetMappingPB generateSortedCallSetMap(
            final Map<String, FeatureReader<VariantContext>> sampleToReaderMap,
            final boolean useSamplesInOrderProvided,
            final boolean validateSampleMap,
            final long lbRowIdx) {
        if (!validateSampleMap) return generateSortedCallSetMap(new ArrayList<>(sampleToReaderMap.keySet()),
                useSamplesInOrderProvided, lbRowIdx);

        List<String> listOfSampleNames = new ArrayList<>(sampleToReaderMap.size());
        for (Map.Entry<String, FeatureReader<VariantContext>> mapObject : sampleToReaderMap.entrySet()) {
            String sampleName = mapObject.getKey();
            FeatureReader<VariantContext> featureReader = mapObject.getValue();

            if (featureReader == null)
                throw new IllegalArgumentException("Null FeatureReader found for sample: " + sampleName);

            VCFHeader header = (VCFHeader) featureReader.getHeader();
            List<String> sampleNamesInHeader = header.getSampleNamesInOrder();
            if (sampleNamesInHeader.size() > 1) {
                StringBuilder messageBuilder = new StringBuilder("Multiple samples ");
                for (String name : sampleNamesInHeader)
                    messageBuilder.append(name).append(" ");
                messageBuilder.append(" appear in header for ").append(sampleName);
                throw new IllegalArgumentException(messageBuilder.toString());
            } else {
                if (!sampleName.equals(sampleNamesInHeader.get(0))) {
                    System.err.println("Sample " + header + " does not match with " + sampleNamesInHeader.get(0) +
                            " in header");
                }
            }
            listOfSampleNames.add(sampleName);
        }
        return generateSortedCallSetMap(listOfSampleNames, useSamplesInOrderProvided, lbRowIdx);
    }
}
