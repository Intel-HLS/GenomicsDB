package com.intel.genomicsdb;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

import java.util.Iterator;

/**
 * Utility class wrapping a stream and a VariantContextWriter for a given stream
 * Each GenomicsDB import stream consists of a buffer stream and a writer object
 * If the caller provides an iterator, then mCurrentVC points to the VariantContext
 * object to be written, if any.
 */
class GenomicsDBImporterStreamWrapper
{
  VariantContextWriter mVCWriter = null;
  SilentByteBufferStream mStream = null;
  private Iterator<VariantContext> mIterator = null;
  private VariantContext mCurrentVC = null;

  /**
   * Constructor
   * @param vcfHeader VCF header for the stream
   * @param bufferCapacity Capacity of the stream buffer in bytes
   * @param streamType BCF_STREAM or VCF_STREAM
   * @param vcIterator iterator over VariantContext objects,
   *                   can be null if the caller is managing
   *                   the buffer explicitly
   */
  public GenomicsDBImporterStreamWrapper(
    final VCFHeader vcfHeader,
    final long bufferCapacity,
    final VariantContextWriterBuilder.OutputType streamType,
    Iterator<VariantContext> vcIterator) throws GenomicsDBException
  {
    mIterator = vcIterator;
    if(vcIterator != null && vcIterator.hasNext())
      mCurrentVC = vcIterator.next();
    boolean headerWritten = false;
    long currentCapacity = bufferCapacity;
    //Must ensure that the header gets written into the buffer stream
    //Why this big outer loop? VCFWriter/BCFWriter seems to store some state which makes
    //calling writeHeader() multiple times impossible
    //Hence, create new objects in every iteration of the loop
    //Since this function is called only once per stream, not really
    //a performance concern
    while(!headerWritten)
    {
      mStream = new SilentByteBufferStream(currentCapacity);
      switch(streamType)
      {
        case BCF_STREAM:
          mVCWriter = new VariantContextWriterBuilder()
            .setOutputBCFStream(mStream)
            .unsetOption(Options.INDEX_ON_THE_FLY).build();
          break;
        case VCF_STREAM:
          mVCWriter = new VariantContextWriterBuilder()
            .setOutputVCFStream(mStream)
            .unsetOption(Options.INDEX_ON_THE_FLY).build();
          break;
        default:
          throw new GenomicsDBException("Unknown stream type "+streamType.toString());
      }
      //Why clone the header?
      //The writer modifies the VCFHeader object passed to writeHeader() - however,
      //we might need to call writeHeader multiple times if the underlying buffer
      //in mStream is too small. Hence, always pass a clone of the original,
      //unmodified header in each call of writeHeader
      mVCWriter.writeHeader(new VCFHeader(vcfHeader));
      if(mStream.overflow())
        currentCapacity = 2*currentCapacity+1;
      else
        headerWritten = true;
    }
  }

  /**
   * Returns true if a non-null Iterator over VariantContext
   * objects was provided for this stream
   * @return True if iterator was provided, False otherwise
   */
  boolean hasIterator()
  {
    return (mIterator != null);
  }

  /**
   * Returns the next VariantContext object iff the Iterator over VariantContext objects
   * is non-null and has a next() object,
   * else returns null. Stores the result in mCurrentVC
   * @return the next VariantContext object or null
   */
  public VariantContext next()
  {
    if(mIterator != null && mIterator.hasNext())
      mCurrentVC = mIterator.next();
    else
      mCurrentVC = null;
    return mCurrentVC;
  }

  /**
   * Returns mCurrentVC - could be null if mIterator is null or !mIterator.hasNext()
   * @return VariantContext object to be written
   */
  VariantContext getCurrentVC()
  {
    return mCurrentVC;
  }
}