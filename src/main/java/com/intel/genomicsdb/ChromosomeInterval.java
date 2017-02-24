package com.intel.genomicsdb;

/**
 * Utility class to represent a chromosome interval
 * Contains 3 members - chr name, start, end (1-based)
 */
class ChromosomeInterval
{
  String mChromosomeName = null;
  long mBegin = -1;
  long mEnd = -1;

  ChromosomeInterval(final String name, final long begin, final long end)
  {
    mChromosomeName = name;
    mBegin = begin;
    mEnd = end;
  }
}
