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

package com.intel.genomicsdb.importer.model;

import htsjdk.samtools.util.Locatable;

/**
 * Utility class to represent a chromosome interval
 * Contains 3 members - chr name, start, end (1-based)
 */
public class ChromosomeInterval implements Locatable {

  private String mChromosomeName = null;
  private long mBegin = -1;
  private long mEnd = -1;

  public ChromosomeInterval(final String name, final long begin, final long end) {
    mChromosomeName = name;
    mBegin = begin;
    mEnd = end;
  }

  @Override
  public String getContig() {
    return mChromosomeName;
  }

  @Override
  public int getStart() {
    return (int) mBegin;
  }

  @Override
  public int getEnd() {
    return (int) mEnd;
  }
}
