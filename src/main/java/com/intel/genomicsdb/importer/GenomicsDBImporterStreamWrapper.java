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

import com.intel.genomicsdb.exception.GenomicsDBException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

import java.util.Iterator;

/**
 * Utility class wrapping a stream and a VariantContextWriter for a given stream
 * Each GenomicsDB import stream consists of a buffer stream and a writer object
 * If the caller provides an iterator, then currentVC points to the VariantContext
 * object to be written, if any.
 */
class GenomicsDBImporterStreamWrapper {
    private Iterator<VariantContext> iterator;
    private VariantContext currentVC = null;
    private VariantContextWriter vcWriter = null;
    private SilentByteBufferStream stream = null;

    /**
     * Constructor
     *
     * @param vcfHeader      VCF header for the stream
     * @param bufferCapacity Capacity of the stream buffer in bytes
     * @param streamType     BCF_STREAM or VCF_STREAM
     * @param vcIterator     iterator over VariantContext objects,
     *                       can be null if the caller is managing
     *                       the buffer explicitly
     */
    public GenomicsDBImporterStreamWrapper(final VCFHeader vcfHeader, final long bufferCapacity,
                                           final VariantContextWriterBuilder.OutputType streamType,
                                           Iterator<VariantContext> vcIterator) throws GenomicsDBException {
        iterator = vcIterator;
        if (vcIterator != null && vcIterator.hasNext()) currentVC = vcIterator.next();
        boolean headerWritten = false;
        long currentCapacity = bufferCapacity;
        //Must ensure that the header gets written into the buffer stream
        //Why this big outer loop? VCFWriter/BCFWriter seems to store some state which makes
        //calling writeHeader() multiple times impossible
        //Hence, create new objects in every iteration of the loop
        //Since this function is called only once per stream, not really
        //a performance concern
        while (!headerWritten) {
            this.stream = new SilentByteBufferStream(currentCapacity);
            switch (streamType) {
                case BCF_STREAM:
                    this.vcWriter = new VariantContextWriterBuilder().setOutputBCFStream(getStream())
                            .unsetOption(Options.INDEX_ON_THE_FLY).build();
                    break;
                case VCF_STREAM:
                    this.vcWriter = new VariantContextWriterBuilder().setOutputVCFStream(getStream())
                            .unsetOption(Options.INDEX_ON_THE_FLY).build();
                    break;
                default:
                    throw new GenomicsDBException("Unknown stream type " + streamType.toString());
            }
            //Why clone the header?
            //The writer modifies the VCFHeader object passed to writeHeader() - however,
            //we might need to call writeHeader multiple times if the underlying buffer
            //in stream is too small. Hence, always pass a clone of the original,
            //unmodified header in each call of writeHeader
            getVcWriter().writeHeader(new VCFHeader(vcfHeader));
            if (getStream().overflow()) currentCapacity = 2 * currentCapacity + 1;
            else headerWritten = true;
        }
    }

    /**
     * Returns true if a non-null Iterator over VariantContext
     * objects was provided for this stream
     *
     * @return True if iterator was provided, False otherwise
     */
    boolean hasIterator() { return iterator != null; }

    /**
     * Returns the next VariantContext object iff the Iterator over VariantContext objects
     * is non-null and has a next() object,
     * else returns null. Stores the result in currentVC
     *
     * @return the next VariantContext object or null
     */
    public VariantContext next() {
        if (iterator != null && iterator.hasNext()) currentVC = iterator.next();
        else currentVC = null;
        return currentVC;
    }

    /**
     * Returns currentVC - could be null if iterator is null or !iterator.hasNext()
     *
     * @return VariantContext object to be written
     */
    VariantContext getCurrentVC() {
        return currentVC;
    }

    public VariantContextWriter getVcWriter() {
        return vcWriter;
    }

    public SilentByteBufferStream getStream() {
        return stream;
    }
}
