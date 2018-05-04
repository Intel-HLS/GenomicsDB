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

import java.io.IOException;
import java.io.OutputStream;

import static com.intel.genomicsdb.importer.Constants.DEFAULT_BUFFER_CAPACITY;

/**
 * Buffer stream implementation - it's silent in the sense that when the buffer is full,
 * it doesn't raise an exception but just marks a a flag as full. It's up to the caller
 * to check the flag and retry later
 * Why? Most likely, it's faster to check a flag rather than throw and catch an exception
 */
class SilentByteBufferStream extends OutputStream {
    private byte buffer[] = null;
    private long numValidBytes = 0;
    private long marker = 0;
    private boolean overflow = false;

    /**
     * Constructor - uses default value of buffer capacity (20KiB)
     */
    public SilentByteBufferStream() {
        buffer = new byte[(int) DEFAULT_BUFFER_CAPACITY];
    }

    /**
     * Constructor - uses specified buffer capacity
     *
     * @param capacity size of buffer in bytes
     */
    public SilentByteBufferStream(final long capacity) {
        buffer = new byte[(int) capacity];
    }

    @Override
    public void close() throws IOException  //does nothing
    {
    }

    @Override
    public void flush() throws IOException  //does nothing
    {
    }

    @Override
    public void write(byte[] b, int off, int len) throws IOException {
        if (overflow) return;
        if (len + numValidBytes > buffer.length) overflow = true;
        else {
            System.arraycopy(b, off, buffer, (int) numValidBytes, len);
            numValidBytes += len;
        }
    }

    @Override
    public void write(byte[] b) throws IOException { write(b, 0, b.length); }

    @Override
    public void write(int b) throws IOException {
        if (overflow) return;
        if (numValidBytes + 1 > buffer.length) overflow = true;
        else {
            buffer[(int) numValidBytes] = (byte) b;
            ++numValidBytes;
        }
    }

    /**
     * Returns buffer capacity in bytes
     *
     * @return buffer capacity in bytes
     */
    public int size() {
        return buffer.length;
    }

    /**
     * Resizes buffer to new size - data is retained
     *
     * @param newSize new capacity of the buffer
     */
    public void resize(final long newSize) {
        byte tmp[] = new byte[(int) newSize];
        System.arraycopy(buffer, 0, tmp, 0, buffer.length);
        buffer = tmp; //hopefully Java GC does its job
    }

    /**
     * Returns if the buffer has overflowed
     *
     * @return true if the buffer has overflowed
     */
    public boolean overflow() {
        return overflow;
    }

    /**
     * Set overflow value
     *
     * @param value overflow value
     */
    public void setOverflow(final boolean value) {
        overflow = value;
    }

    /**
     * Get number of valid bytes
     *
     * @return number of valid bytes
     */
    public long getNumValidBytes() {
        return numValidBytes;
    }

    /**
     * Set number of valid bytes
     *
     * @param value number of valid bytes
     */
    public void setNumValidBytes(final long value) {
        numValidBytes = value;
    }

    /**
     * Get marker value
     *
     * @return marker value
     */
    public long getMarker() {
        return marker;
    }

    /**
     * Caller code can use this function to mark a certain point in the buffer
     * This is generally used to mark the position in the buffer after the last
     * complete VariantContext object written
     *
     * @param value set marker value
     */
    public void setMarker(final long value) {
        marker = value;
    }

    /**
     * Get byte buffer for this stream
     *
     * @return byte buffer for this stream
     */
    public byte[] getBuffer() {
        return buffer;
    }
}
