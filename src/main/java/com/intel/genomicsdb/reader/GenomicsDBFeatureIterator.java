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

import com.intel.genomicsdb.model.Coordinates;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.variant.bcf2.BCF2Codec;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalInt;
import java.util.stream.Collectors;

/**
 * Iterator over {@link htsjdk.variant.variantcontext.VariantContext} objects.
 * Uses {@link GenomicsDBQueryStream} to obtain combined gVCF records
 * (as BCF2) from TileDB/GenomicsDB
 */
public class GenomicsDBFeatureIterator<T extends Feature, SOURCE> implements CloseableTribbleIterator<T> {
    private FeatureCodecHeader featureCodecHeader;
    private FeatureCodec<T, SOURCE> codec;
    private List<SOURCE> sources;
    private GenomicsDBTimer timer;
    private boolean closedBefore;
    private SOURCE currentSource;

    /**
     * Constructor
     *
     * @param loaderJSONFile     GenomicsDB loader JSON configuration file
     * @param queryJSONFiles     GenomicsDB query JSON configuration file list
     * @param featureCodecHeader htsjdk Feature codec header
     * @param codec              FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *                           and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @param intervalPerArray   Optional map of relation between qjf and a contig interval
     * @throws IOException when data cannot be read from the stream
     */
    GenomicsDBFeatureIterator(final String loaderJSONFile, final List<String> queryJSONFiles,
                              final FeatureCodecHeader featureCodecHeader, final FeatureCodec<T, SOURCE> codec,
                              final Optional<Map<String, Coordinates.ContigInterval>> intervalPerArray)
            throws IOException {
        this(loaderJSONFile, queryJSONFiles, featureCodecHeader, codec, "", OptionalInt.empty(),
                OptionalInt.empty(), intervalPerArray);
    }

    /**
     * Constructor
     *
     * @param loaderJSONFile     GenomicsDB loader JSON configuration file
     * @param queryJSONFiles     GenomicsDB query JSON configuration file list
     * @param featureCodecHeader htsjdk Feature codec header
     * @param codec              FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *                           and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @param chr                contig name
     * @param start              start position (1-based)
     * @param end                end position, inclusive (1-based)
     * @param intervalPerArray   Optional map of relation between qjf and a contig interval
     * @throws IOException when data cannot be read from the stream
     */
    GenomicsDBFeatureIterator(final String loaderJSONFile, final List<String> queryJSONFiles,
                              final FeatureCodecHeader featureCodecHeader, final FeatureCodec<T, SOURCE> codec,
                              final String chr, final OptionalInt start, final OptionalInt end,
                              final Optional<Map<String, Coordinates.ContigInterval>> intervalPerArray) throws IOException {
        this.featureCodecHeader = featureCodecHeader;
        this.codec = codec;
        boolean readAsBCF = this.codec instanceof BCF2Codec;
        boolean areIntervalPerArraySpecified = intervalPerArray.isPresent() && intervalPerArray.get().size() > 0;
        this.sources = queryJSONFiles.stream().map(qjf -> {
            GenomicsDBQueryStream genomicsDBQueryStream;
            if (areIntervalPerArraySpecified && start.isPresent() && end.isPresent()) {
                Coordinates.ContigInterval interval = intervalPerArray.get().get(qjf);
                genomicsDBQueryStream = new GenomicsDBQueryStream(loaderJSONFile, qjf, interval.getContig(),
                        Math.max(start.getAsInt(), (int) interval.getBegin()),
                        Math.min(end.getAsInt(), (int) interval.getEnd()), readAsBCF);
            } else if (areIntervalPerArraySpecified && !start.isPresent() && !end.isPresent()) {
                Coordinates.ContigInterval interval = intervalPerArray.get().get(qjf);
                genomicsDBQueryStream = new GenomicsDBQueryStream(loaderJSONFile, qjf, interval.getContig(),
                        (int) interval.getBegin(), (int) interval.getEnd(), readAsBCF);
            } else {
                genomicsDBQueryStream = new GenomicsDBQueryStream(loaderJSONFile, qjf, chr, 0, 0, readAsBCF);
            }
            if (readAsBCF) { //BCF2 codec provides size of header
                try {
                    genomicsDBQueryStream.skip(this.featureCodecHeader.getHeaderEnd());
                } catch (IOException ex) {
                    throw new RuntimeException(ex);
                }
            }
            return genomicsDBQueryStream;
        }).map(gqs -> this.codec.makeSourceFromStream(gqs)).collect(Collectors.toList());
        if (sources.isEmpty()) throw new IllegalStateException("There are no sources based on those query parameters");
        if (!readAsBCF) { //VCF Codec must parse out header again since getHeaderEnd() returns 0
          for(SOURCE currSource : this.sources)
            this.codec.readHeader(currSource); //no need to store header anywhere
        }
        this.currentSource = this.sources.get(0);
        this.timer = new GenomicsDBTimer();
        this.closedBefore = false;
    }

    @Override
    public boolean hasNext() {
        int index = 0;
        //While loop since the next source might not return any data, but subsequent sources might
        while(this.codec.isDone(this.currentSource) && index < this.sources.size())
          index = setNextSourceAsCurrent();
        boolean isDone = (this.codec.isDone(this.currentSource));
        if (isDone) close();
        return !isDone;
    }

    @Override
    public T next() {
        try {
            this.timer.start();
            if (this.codec.isDone(this.currentSource)) throw new RuntimeException("No more valid records exist");
            T nextObj = this.codec.decode(this.currentSource);
            this.timer.stop();
            return nextObj;
        } catch (IOException e) {
            throw new RuntimeException("Unknown codec exception");
        }
    }

    @Override
    public void close() {
        if (!this.closedBefore) {
            this.timer.print("GenomicsDB iterator next() timer", System.err);
            this.codec.close(this.currentSource);
            this.closedBefore = true;
        }
    }

    @Override
    public GenomicsDBFeatureIterator iterator() {
        return this;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported in Iterators");
    }

    private int setNextSourceAsCurrent() {
        this.codec.close(this.currentSource);
        int index = this.sources.indexOf(this.currentSource);
        if (index < this.sources.size() - 1) this.currentSource = this.sources.get(index + 1);
        return index + 1;
    }
}
