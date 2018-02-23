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

import com.intel.genomicsdb.GenomicsDBQueryStream;
import com.intel.genomicsdb.reader.model.QueryParams;
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
import java.util.stream.Stream;

import static java.util.stream.Collectors.toList;

/**
 * A reader for GenomicsDB that implements {@link htsjdk.tribble.FeatureReader}
 * Currently, the reader only return {@link htsjdk.variant.variantcontext.VariantContext}
 */
public class GenomicsDBFeatureReader<T extends Feature, SOURCE> implements FeatureReader<T> {
    private String loaderJSONFile;
    private QueryParams queryParams;
    private FeatureCodec<T, SOURCE> codec;
    private FeatureCodecHeader featureCodecHeader;
    private ArrayList<String> sequenceNames;

    /**
     * Constructor
     *
     * @param loaderJSONFile  GenomicsDB loader JSON configuration file
     * @param tiledbWorkspace TileDB workspace path
     * @param referenceGenome Path to reference genome (fasta file)
     * @param codec           FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *                        and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @throws IOException when data cannot be read from the stream
     */
    public GenomicsDBFeatureReader(final String loaderJSONFile, final String tiledbWorkspace,
                                   final String referenceGenome, final FeatureCodec<T, SOURCE> codec) throws IOException {
        this(loaderJSONFile, tiledbWorkspace, referenceGenome, null, codec);
    }

    /**
     * Constructor
     *
     * @param loaderJSONFile            GenomicsDB loader JSON configuration file
     * @param tiledbWorkspace           TileDB workspace path
     * @param referenceGenome           Path to reference genome (fasta file)
     * @param templateVCFHeaderFilename Template VCF header to be used for
     *                                  the combined gVCF records
     * @param codec                     FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *                                  and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @throws IOException when data cannot be read from the stream
     */
    public GenomicsDBFeatureReader(final String loaderJSONFile, final String tiledbWorkspace,
                                   final String referenceGenome, final String templateVCFHeaderFilename,
                                   final FeatureCodec<T, SOURCE> codec) throws IOException {
        initialize(loaderJSONFile, new QueryParams(tiledbWorkspace, referenceGenome, templateVCFHeaderFilename), codec);
    }

    /**
     * Constructor
     *
     * @param vidMappingFile            GenomicsDB vid mapping JSON configuration file
     * @param callsetMappingFile        GenomicsDB callset mapping JSON configuration file
     * @param tiledbWorkspace           TileDB workspace path
     * @param referenceGenome           Path to reference genome (fasta file)
     * @param templateVCFHeaderFilename Template VCF header to be used for
     *                                  the combined gVCF records
     * @param codec                     FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *                                  and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @throws IOException when data cannot be read from the stream
     */
    public GenomicsDBFeatureReader(final String vidMappingFile, final String callsetMappingFile,
                                   final String tiledbWorkspace, final String referenceGenome,
                                   final String templateVCFHeaderFilename, final FeatureCodec<T, SOURCE> codec) throws IOException {
        this(vidMappingFile, callsetMappingFile, tiledbWorkspace, referenceGenome, templateVCFHeaderFilename, codec, false);
    }

    /**
     * Constructor
     *
     * @param vidMappingFile            GenomicsDB vid mapping JSON configuration file
     * @param callsetMappingFile        GenomicsDB callset mapping JSON configuration file
     * @param tiledbWorkspace           TileDB workspace path
     * @param referenceGenome           Path to reference genome (fasta file)
     * @param templateVCFHeaderFilename Template VCF header to be used for
     *                                  the combined gVCF records
     * @param codec                     FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *                                  and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @param produceGT                 By default, GenomicsDB ignores the GT field to match the output of GATK
     *                                  CombineGVCFs. Enabling this flag will produce the GT fields
     * @throws IOException when data cannot be read from the stream
     */
    public GenomicsDBFeatureReader(final String vidMappingFile, final String callsetMappingFile,
                                   final String tiledbWorkspace, final String referenceGenome,
                                   final String templateVCFHeaderFilename, final FeatureCodec<T, SOURCE> codec,
                                   final boolean produceGT) throws IOException {
        initialize("", new QueryParams(tiledbWorkspace, vidMappingFile, callsetMappingFile, referenceGenome,
                templateVCFHeaderFilename, produceGT), codec);
    }

    /**
     * Initialization function that's used by all constructors
     *
     * @param loaderJSONFile GenomicsDB loader JSON configuration file
     * @param queryParams    GenomicsDB query parameters
     * @param codec          FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *                       and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @throws IOException when data cannot be read from the stream
     */
    public void initialize(final String loaderJSONFile, final QueryParams queryParams,
                           final FeatureCodec<T, SOURCE> codec) throws IOException {
        this.loaderJSONFile = loaderJSONFile;
        this.queryParams = queryParams;
        this.codec = codec;
        String[] chromosomeIntervalArrays = getArrayListFromWorkspace(new File(queryParams.getWorkspace()));
        if (chromosomeIntervalArrays.length < 1)
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
        List<String> chromosomeIntervalArraysPaths = resolveChromosomeArrayFolderList();
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
        List<String> chromosomeIntervalArraysPaths = resolveChromosomeArrayFolderList(chr, start, end);
        return new GenomicsDBFeatureIterator(this.loaderJSONFile, chromosomeIntervalArraysPaths, this.featureCodecHeader,
                this.codec, chr, start, end);
    }

    private List<String> resolveChromosomeArrayFolderList() {
        String[] chromosomeIntervalArraysNames = getArrayListFromWorkspace(new File(queryParams.getWorkspace()));
        return createArrayFolderListFromArrayStream(Arrays.stream(chromosomeIntervalArraysNames));
    }

    private List<String> resolveChromosomeArrayFolderList(final String chromosome, final int intervalStart, final int intervalEnd) {
        String[] chromosomeIntervalArraysNames = getArrayListFromWorkspace(new File(queryParams.getWorkspace()));
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
            QueryParams fullQueryParams = new QueryParams(this.queryParams);
            fullQueryParams.setArrayName(name);
            try {
                return createTempQueryJsonFile(fullQueryParams.getArrayName(), fullQueryParams).getAbsolutePath();
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }).collect(toList());
    }

    private void generateHeadersForQuery(final String randomExistingArrayName) throws IOException {
        QueryParams fullQueryParams = new QueryParams(this.queryParams);
        fullQueryParams.setArrayName(randomExistingArrayName);
        File queryJSONFile = createTempQueryJsonFile(randomExistingArrayName, fullQueryParams);
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
    private File createTempQueryJsonFile(final String arrayName, final QueryParams fullQueryParams) throws IOException {
        File tmpQueryJSONFile = File.createTempFile(String.format("queryJSON_%s", arrayName), ".json");
        tmpQueryJSONFile.deleteOnExit();
        FileWriter fptr = new FileWriter(tmpQueryJSONFile);
        fptr.write(fullQueryParams.toJsonString());
        fptr.close();
        return tmpQueryJSONFile;
    }

    private String[] getArrayListFromWorkspace(final File workspace) {
        return workspace.list((current, name) -> new File(current, name).isDirectory());
    }
}