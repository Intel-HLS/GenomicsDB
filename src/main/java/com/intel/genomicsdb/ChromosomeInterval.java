package com.intel.genomicsdb;

import htsjdk.samtools.util.Locatable;

/**
 * Utility class to represent a chromosome interval
 * Contains 3 members - chr name, start, end (1-based)
 */
public class ChromosomeInterval {

  String mChromosomeName = null;
  long mBegin = -1;
  long mEnd = -1;

  public ChromosomeInterval(final String name, final long begin, final long end)
  {
    mChromosomeName = name;
    mBegin = begin;
    mEnd = end;
  }
}
