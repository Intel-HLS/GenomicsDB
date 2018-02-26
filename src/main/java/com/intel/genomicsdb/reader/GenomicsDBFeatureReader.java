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
import com.intel.genomicsdb.GenomicsDBExportConfiguration;
import com.intel.genomicsdb.GenomicsDBQueryStream;
import htsjdk.tribble.*;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

import static java.util.stream.Collectors.toList;

/**
 * A reader for GenomicsDB that implements {@link htsjdk.tribble.FeatureReader}
 * Currently, the reader only return {@link htsjdk.variant.variantcontext.VariantContext}
 */
public class GenomicsDBFeatureReader<T extends Feature, SOURCE> implements FeatureReader<T> {
    private String loaderJSONFile;
    private GenomicsDBExportConfiguration.ExportConfiguration exportConfiguration;
    private FeatureCodec<T, SOURCE> codec;
    private FeatureCodecHeader featureCodecHeader;
    private ArrayList<String> sequenceNames;

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
        String[] chromosomeIntervalArrays = this.exportConfiguration.hasArray() ? new String[]{
                exportConfiguration.getArray()
        } : getArrayListFromWorkspace(new File(exportConfiguration.getWorkspace()));
        if (chromosomeIntervalArrays == null || chromosomeIntervalArrays.length < 1)
            throw new IllegalStateException("There is no genome data stored in the database");
        generateHeadersForQuery(chromosomeIntervalArrays[0]);
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
        List<String> chromosomeIntervalArraysPaths = this.exportConfiguration.hasArray() ? createArrayFolderListFromArrayStream(
                new ArrayList<String>() {{
                    add(exportConfiguration.getArray());
                }}.stream()) : resolveChromosomeArrayFolderList();
        return new GenomicsDBFeatureIterator(this.loaderJSONFile, chromosomeIntervalArraysPaths,
                this.featureCodecHeader, this.codec);
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
        List<String> chromosomeIntervalArraysPaths = this.exportConfiguration.hasArray() ? createArrayFolderListFromArrayStream(
                new ArrayList<String>() {{
                    add(exportConfiguration.getArray());
                }}.stream()) : resolveChromosomeArrayFolderList(chr, start, end);
        return new GenomicsDBFeatureIterator(this.loaderJSONFile, chromosomeIntervalArraysPaths, this.featureCodecHeader,
                this.codec, chr, start, end);
    }

    private List<String> resolveChromosomeArrayFolderList() {
        String[] chromosomeIntervalArraysNames = getArrayListFromWorkspace(new File(exportConfiguration.getWorkspace()));
        return createArrayFolderListFromArrayStream(Arrays.stream(chromosomeIntervalArraysNames));
    }

    private List<String> resolveChromosomeArrayFolderList(final String chromosome, final int intervalStart, final int intervalEnd) {
        String[] chromosomeIntervalArraysNames = getArrayListFromWorkspace(new File(exportConfiguration.getWorkspace()));
        Stream<String> stream = Arrays.stream(chromosomeIntervalArraysNames).filter(name -> {
            String[] ref = name.split("_");
            if (ref.length != 3) throw new RuntimeException("There is a wrong array folder name in the workspace. " +
                    "Array folder name should be {chromosome}_{intervalStart}_{intervalEnd}");
            return chromosome.equals(ref[0]) && ((intervalStart >= Integer.parseInt(ref[1]) && intervalStart <= Integer.parseInt(ref[2])) ||
                    (intervalEnd >= Integer.parseInt(ref[1]) && intervalEnd <= Integer.parseInt(ref[2])));
        });
        return createArrayFolderListFromArrayStream(stream);
    }

    private List<String> createArrayFolderListFromArrayStream(Stream<String> stream) {
        return stream.map(name -> {
            GenomicsDBExportConfiguration.ExportConfiguration fullExportConfiguration =
                    GenomicsDBExportConfiguration.ExportConfiguration.newBuilder(this.exportConfiguration).setArray(name).build();
            try {
                return createTempQueryJsonFile(fullExportConfiguration.getArray(), fullExportConfiguration).getAbsolutePath();
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }).collect(toList());
    }

    private void generateHeadersForQuery(final String randomExistingArrayName) throws IOException {
        GenomicsDBExportConfiguration.ExportConfiguration fullExportConfiguration =
                GenomicsDBExportConfiguration.ExportConfiguration.newBuilder(this.exportConfiguration)
                        .setArray(randomExistingArrayName).build();
        File queryJSONFile = createTempQueryJsonFile(randomExistingArrayName, fullExportConfiguration);
        GenomicsDBQueryStream gdbStream = new GenomicsDBQueryStream(this.loaderJSONFile, queryJSONFile.getAbsolutePath(),
                this.codec instanceof BCF2Codec, true);
        SOURCE source = this.codec.makeSourceFromStream(gdbStream);
        this.featureCodecHeader = this.codec.readHeader(source);
        //Store sequence names
        VCFHeader vcfHeader = (VCFHeader) (this.featureCodecHeader.getHeaderValue());
        this.sequenceNames = new ArrayList<>(vcfHeader.getContigLines().size());
        for (final VCFContigHeaderLine contigHeaderLine : vcfHeader.getContigLines())
            this.sequenceNames.add(contigHeaderLine.getID());
        gdbStream.close();
        queryJSONFile.delete();
    }

    // TODO: remove this once protobuf classes are created
    private File createTempQueryJsonFile(final String arrayName,
                                         final GenomicsDBExportConfiguration.ExportConfiguration finalExportConfiguration) throws IOException {
        File tmpQueryJSONFile = File.createTempFile(String.format("queryJSON_%s", arrayName), ".json");
        tmpQueryJSONFile.deleteOnExit();
        FileWriter fptr = new FileWriter(tmpQueryJSONFile);
        String jsonString = JsonFormat.printToString(finalExportConfiguration);
        fptr.write(jsonString);
        fptr.close();
        return tmpQueryJSONFile;
    }

    private String[] getArrayListFromWorkspace(final File workspace) {
        return workspace.list((current, name) -> new File(current, name).isDirectory());
    }
}