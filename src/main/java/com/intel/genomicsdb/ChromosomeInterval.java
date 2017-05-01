package com.intel.genomicsdb;

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
