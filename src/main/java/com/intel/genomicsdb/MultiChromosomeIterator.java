package com.intel.genomicsdb;

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
 * @param <SOURCE> LineIterator for VCFs, PositionalBufferedStream for BCFs
 */
class MultiChromosomeIterator<SOURCE> implements Iterator<VariantContext>
{
  private ArrayList<ChromosomeInterval> mChromosomeIntervals = null;
  private AbstractFeatureReader<VariantContext, SOURCE> mReader = null;
  private int mIdxInIntervalList = 0;
  private CloseableTribbleIterator<VariantContext> mIterator = null;

  /**
   * Constructor
   * @param reader AbstractFeatureReader over VariantContext objects -
   *               SOURCE can vary - BCF v/s VCF for example
   * @param chromosomeIntervals chromosome intervals over which to iterate
   * @throws IOException when the reader's query method throws an IOException
   */
  MultiChromosomeIterator(AbstractFeatureReader<VariantContext, SOURCE> reader,
                          final List<ChromosomeInterval> chromosomeIntervals)
    throws IOException
  {

    mReader = reader;
    VCFHeader mHeader = (VCFHeader) (reader.getHeader());
    mChromosomeIntervals = new ArrayList<>();
    //Only add intervals whose chromosomes are present in the VCF header
    final SAMSequenceDictionary contigDictionary = mHeader.getSequenceDictionary();
    for(ChromosomeInterval currInterval : chromosomeIntervals)
      if(contigDictionary.getSequenceIndex(currInterval.getContig()) != -1)
        mChromosomeIntervals.add(currInterval);
    if(mChromosomeIntervals.size() > 0)
    {
      ChromosomeInterval currInterval = mChromosomeIntervals.get(0);
      mIterator = mReader.query(currInterval.getContig(),
        (int)currInterval.getStart(),
        (int)currInterval.getEnd());
    }
  }

  @Override
  public boolean hasNext() {
    return mIterator != null && mIterator.hasNext();
  }

  @Override
  public VariantContext next() throws NoSuchElementException
  {
    try
    {
      if(mIterator == null)
        throw new NoSuchElementException("next() called for iterator with no more elements");
      VariantContext returnValue = mIterator.next();
      //within the same chromosome
      if(mIterator.hasNext())
        return returnValue;
      //move to next chromosome and iterate
      //It's possible that the reader has no record for the next contig, but could have records
      //for subsequent contigs
      for(mIdxInIntervalList=mIdxInIntervalList+1;mIdxInIntervalList<mChromosomeIntervals.size();
          ++mIdxInIntervalList)
      {
        ChromosomeInterval currInterval = mChromosomeIntervals.get(mIdxInIntervalList);
        mIterator = mReader.query(currInterval.getContig(),
          (int)currInterval.getStart(),
          (int)currInterval.getEnd());
        if(mIterator.hasNext())
          return returnValue;
      }
      mIterator = null;
      return returnValue;
    }
    catch(IOException e)
    {
      throw new NoSuchElementException("Caught IOException: "+e.getMessage());
    }
  }
} // End of class MultiChromosomeIterator
