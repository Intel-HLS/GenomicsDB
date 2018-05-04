/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

package com.intel.genomicsdb.importer;

import com.intel.genomicsdb.importer.model.ChromosomeInterval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Given an AbstractFeatureReader over VariantContext objects
 * and a list of ChromosomeInterval objects, this class
 * is an Iterator over VariantContext for all the chromosome intervals in the
 * list
 *
 * @param <SOURCE> LineIterator for VCFs, PositionalBufferedStream for BCFs
 */
class MultiChromosomeIterator<SOURCE> implements Iterator<VariantContext> {
    private ArrayList<ChromosomeInterval> chromosomeIntervals;
    private AbstractFeatureReader<VariantContext, SOURCE> reader;
    private int idxInIntervalList = 0;
    private CloseableTribbleIterator<VariantContext> iterator = null;

    /**
     * Constructor
     *
     * @param reader              AbstractFeatureReader over VariantContext objects -
     *                            SOURCE can vary - BCF v/s VCF for example
     * @param chromosomeIntervals chromosome intervals over which to iterate
     * @throws IOException when the reader's query method throws an IOException
     */
    MultiChromosomeIterator(AbstractFeatureReader<VariantContext, SOURCE> reader, final List<ChromosomeInterval> chromosomeIntervals)
            throws IOException {
        this.reader = reader;
        this.chromosomeIntervals = new ArrayList<>();
        VCFHeader header = (VCFHeader) (reader.getHeader());
        //Only add intervals whose chromosomes are present in the VCF header
        final SAMSequenceDictionary contigDictionary = header.getSequenceDictionary();
        for (ChromosomeInterval currInterval : chromosomeIntervals)
            if (contigDictionary.getSequenceIndex(currInterval.getContig()) != -1) this.chromosomeIntervals.add(currInterval);
        if (this.chromosomeIntervals.size() > 0) {
            ChromosomeInterval currInterval = this.chromosomeIntervals.get(0);
            iterator = this.reader.query(currInterval.getContig(), currInterval.getStart(), currInterval.getEnd());
        }
    }

    @Override
    public boolean hasNext() {
        return iterator != null && iterator.hasNext();
    }

    @Override
    public VariantContext next() throws NoSuchElementException {
        try {
            if (iterator == null) throw new NoSuchElementException("next() called for iterator with no more elements");
            VariantContext returnValue = iterator.next();
            //within the same chromosome
            if (iterator.hasNext()) return returnValue;
            //move to next chromosome and iterate
            //It's possible that the reader has no record for the next contig, but could have records
            //for subsequent contigs
            for (idxInIntervalList = idxInIntervalList + 1; idxInIntervalList < chromosomeIntervals.size();
                 ++idxInIntervalList) {
                ChromosomeInterval currInterval = chromosomeIntervals.get(idxInIntervalList);
                iterator = reader.query(currInterval.getContig(), currInterval.getStart(), currInterval.getEnd());
                if (iterator.hasNext()) return returnValue;
            }
            iterator = null;
            return returnValue;
        } catch (IOException e) {
            throw new NoSuchElementException("Caught IOException: " + e.getMessage());
        }
    }
} // End of class MultiChromosomeIterator
