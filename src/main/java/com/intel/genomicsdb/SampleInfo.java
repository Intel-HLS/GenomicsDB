package com.intel.genomicsdb;

/**
 * Utility class that stores row index and globally unique name for a given sample
 */
public class SampleInfo
{
  String mName = null;
  long mRowIdx = -1;

  public SampleInfo(final String name, final long rowIdx)
  {
    mName = name;
    mRowIdx = rowIdx;
  }
}