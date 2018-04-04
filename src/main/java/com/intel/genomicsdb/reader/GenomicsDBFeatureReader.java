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

package com.intel.genomicsdb.reader;

import com.googlecode.protobuf.format.JsonFormat;
import com.intel.genomicsdb.model.Coordinates;
import com.intel.genomicsdb.model.GenomicsDBExportConfiguration;
import htsjdk.tribble.*;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Stream;

import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toMap;

/**
 * A reader for GenomicsDB that implements {@link htsjdk.tribble.FeatureReader}
 * Currently, the reader only return {@link htsjdk.variant.variantcontext.VariantContext}
 */
public class GenomicsDBFeatureReader<T extends Feature, SOURCE> implements FeatureReader<T> {
    private String loaderJSONFile;
    private String queryJsonFileName = null;
    private GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration;
    private FeatureCodec<T, SOURCE> codec;
    private FeatureCodecHeader featureCodecHeader;
    private List<String> sequenceNames;
    private Map<String, Coordinates.ContigInterval> intervalsPerArray = new HashMap<>();

    /**
     * Constructor
     *
     * @param exportConfiguration query parameters
     * @param codec               FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *                            and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @param loaderJSONFile      GenomicsDB loader JSON configuration file
     * @throws IOException when data cannot be read from the stream
     */
    public GenomicsDBFeatureReader(final GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration,
                                   final FeatureCodec<T, SOURCE> codec,
                                   final Optional<String> loaderJSONFile) throws IOException {
        this.exportConfiguration = exportConfiguration;
        this.codec = codec;
        this.loaderJSONFile = loaderJSONFile.orElse("");
        List<String> chromosomeIntervalArrays = this.exportConfiguration.hasArrayName() ? new ArrayList<String>() {{
            add(exportConfiguration.getArrayName());
        }} : getArrayListFromWorkspace(new File(exportConfiguration.getWorkspace()), Optional.empty());
        if (chromosomeIntervalArrays == null || chromosomeIntervalArrays.size() < 1)
            throw new IllegalStateException("There is no genome data stored in the database");
        generateHeadersForQuery(chromosomeIntervalArrays.get(0));
    }

    /**
     * Constructor
     *
     * @param queryJSONFilename JSON file with query parameters
     * @param codec             FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *                          and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @param loaderJSONFile    GenomicsDB loader JSON configuration file
     * @throws IOException when data cannot be read from the stream
     */
    public GenomicsDBFeatureReader(final String queryJSONFilename,
                                   final FeatureCodec<T, SOURCE> codec,
                                   final Optional<String> loaderJSONFile) throws IOException {
        this.queryJsonFileName = queryJSONFilename;
        this.codec = codec;
        this.loaderJSONFile = loaderJSONFile.orElse("");
        generateHeadersForQueryGivenQueryJSONFile(queryJSONFilename);
    }

    /**
     * Return the VCF header of the combined gVCF stream
     *
     * @return the VCF header of the combined gVCF stream
     */
    public Object getHeader() {
        return this.featureCodecHeader.getHeaderValue();
    }

    /**
     * Return the list of contigs in the combined VCF header
     *
     * @return list of strings of the contig names
     */
    public List<String> getSequenceNames() {
        return this.sequenceNames;
    }

    public void close() throws IOException {
    }

    /**
     * Return an iterator over {@link htsjdk.variant.variantcontext.VariantContext}
     * objects for the specified TileDB array and query configuration
     *
     * @return iterator over {@link htsjdk.variant.variantcontext.VariantContext} objects
     */
    public CloseableTribbleIterator<T> iterator() throws IOException {
        List<String> chromosomeIntervalArraysPaths =
                (this.queryJsonFileName != null && !this.queryJsonFileName.isEmpty())
                        ? Collections.singletonList(this.queryJsonFileName)
                        : this.exportConfiguration.hasArrayName() ? createArrayFolderListFromArrayStream(
                        new ArrayList<String>() {{
                            add(exportConfiguration.getArrayName());
                        }}.stream()) : resolveChromosomeArrayFolderList(Optional.empty());
        return new GenomicsDBFeatureIterator(this.loaderJSONFile, chromosomeIntervalArraysPaths,
                this.featureCodecHeader, this.codec, Optional.of(this.intervalsPerArray));
    }

    /**
     * Return an iterator over {@link htsjdk.variant.variantcontext.VariantContext}
     * objects for the specified TileDB array and queried position
     *
     * @param chr   contig name
     * @param start start position (1-based)
     * @param end   end position, inclusive (1-based)
     * @return iterator over {@link htsjdk.variant.variantcontext.VariantContext} objects
     */
    public CloseableTribbleIterator<T> query(final String chr, final int start, final int end) throws IOException {
        List<String> chromosomeIntervalArraysPaths =
                (this.queryJsonFileName != null && !this.queryJsonFileName.isEmpty())
                        ? Collections.singletonList(this.queryJsonFileName)
                        : this.exportConfiguration.hasArrayName() ? createArrayFolderListFromArrayStream(
                        new ArrayList<String>() {{
                            add(exportConfiguration.getArrayName());
                        }}.stream()) : resolveChromosomeArrayFolderList(Optional.of(
                                Coordinates.ContigInterval.newBuilder().setContig(chr).setBegin(start).setEnd(end).build()
                ));
        return new GenomicsDBFeatureIterator(this.loaderJSONFile, chromosomeIntervalArraysPaths, this.featureCodecHeader,
                this.codec, chr, OptionalInt.of(start), OptionalInt.of(end), Optional.of(this.intervalsPerArray));
    }

    private List<String> resolveChromosomeArrayFolderList(final Optional<Coordinates.ContigInterval> chromosome) {
        List<String> chromosomeIntervalArraysNames = getArrayListFromWorkspace(new File(exportConfiguration.getWorkspace()), chromosome);
        chromosomeIntervalArraysNames.sort(new ChrArrayFolderComparator());
        List<String> qjfs = createArrayFolderListFromArrayStream(chromosomeIntervalArraysNames.stream());
        List<Coordinates.ContigInterval> contigIntervals = chromosomeIntervalArraysNames.stream().map(name -> {
            String[] ref = name.split("#");
            throwExceptionIfArrayFolderRefIsWrong(ref);
            return Coordinates.ContigInterval.newBuilder().setContig(ref[0]).setBegin(Integer.parseInt(ref[1]))
                    .setEnd(Integer.parseInt(ref[2])).build();
        }).collect(toList());
        this.intervalsPerArray = contigIntervals.stream().map(contig -> {
            List<String> qjf = qjfs.stream().filter(file -> file.contains(String.format(
                    "%s#%d#%d", contig.getContig(), contig.getBegin(), contig.getEnd()))).collect(toList());
            if (qjf.size() != 1) {
                throw new IllegalStateException("Multiple results for relation array folder name <-> contig interval values");
            }
            return new AbstractMap.SimpleEntry<>(qjf.get(0), contig);
        }).collect(toMap(AbstractMap.SimpleEntry::getKey, AbstractMap.SimpleEntry::getValue));
        return qjfs;
    }

    private void throwExceptionIfArrayFolderRefIsWrong(String[] ref) {
        if (ref.length != 3) throw new RuntimeException("There is a wrong array folder name in the workspace. " +
                "Array folder name should be {chromosome}#{intervalStart}#{intervalEnd}");
    }

    private List<String> createArrayFolderListFromArrayStream(Stream<String> stream) {
        return stream.map(name -> {
            GenomicsDBExportConfiguration.ExportConfiguration fullExportConfiguration =
                    GenomicsDBExportConfiguration.ExportConfiguration.newBuilder(this.exportConfiguration).setArrayName(name).build();
            try {
                return createTempQueryJsonFile(fullExportConfiguration.getArrayName(), fullExportConfiguration).getAbsolutePath();
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }).collect(toList());
    }

    private void generateHeadersForQuery(final String randomExistingArrayName) throws IOException {
        GenomicsDBExportConfiguration.ExportConfiguration fullExportConfiguration =
                GenomicsDBExportConfiguration.ExportConfiguration.newBuilder(this.exportConfiguration)
                        .setArrayName(randomExistingArrayName).build();
        File queryJSONFile = createTempQueryJsonFile(randomExistingArrayName, fullExportConfiguration);
        generateHeadersForQueryGivenQueryJSONFile(queryJSONFile.getAbsolutePath());
        queryJSONFile.delete();
    }

    private void generateHeadersForQueryGivenQueryJSONFile(final String queryJSONFilename) throws IOException {
        GenomicsDBQueryStream gdbStream = new GenomicsDBQueryStream(this.loaderJSONFile, queryJSONFilename,
                this.codec instanceof BCF2Codec, true);
        SOURCE source = this.codec.makeSourceFromStream(gdbStream);
        this.featureCodecHeader = this.codec.readHeader(source);
        //Store sequence names
        VCFHeader vcfHeader = (VCFHeader) (this.featureCodecHeader.getHeaderValue());
        this.sequenceNames = new ArrayList<>(vcfHeader.getContigLines().size());
        for (final VCFContigHeaderLine contigHeaderLine : vcfHeader.getContigLines())
            this.sequenceNames.add(contigHeaderLine.getID());
        gdbStream.close();
    }

    // TODO: remove this once protobuf classes are created
    private File createTempQueryJsonFile(final String arrayName,
                                         final GenomicsDBExportConfiguration.ExportConfiguration finalExportConfiguration)
            throws IOException {
        File tmpQueryJSONFile = File.createTempFile(String.format("queryJSON_%s", arrayName), ".json");
        tmpQueryJSONFile.deleteOnExit();
        FileWriter fptr = new FileWriter(tmpQueryJSONFile);
        String jsonString = JsonFormat.printToString(finalExportConfiguration);
        fptr.write(jsonString);
        fptr.close();
        return tmpQueryJSONFile;
    }

    private List<String> getArrayListFromWorkspace(final File workspace, Optional<Coordinates.ContigInterval> chromosome) {
        List<String> folders = Arrays.asList(workspace.list((current, name) -> new File(current, name).isDirectory()));
        return chromosome.isPresent() ? folders.stream().filter(name -> {
            String[] ref = name.split("#");
            throwExceptionIfArrayFolderRefIsWrong(ref);
            return chromosome.get().getContig().equals(ref[0]) && (chromosome.get().getBegin() <= Integer.parseInt(ref[2])
                    && chromosome.get().getEnd() >= Integer.parseInt(ref[1]));
        }).collect(toList()) : folders;
    }
}
