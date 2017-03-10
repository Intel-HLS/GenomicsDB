import com.intel.genomicsdb.ChromosomeInterval;
import com.intel.genomicsdb.GenomicsDBException;
import com.intel.genomicsdb.GenomicsDBImporter;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.json.simple.parser.ParseException;

import java.io.IOException;
import java.util.*;

/**
 * Created by kdatta1 on 3/10/17.
 */
public final class ExampleGenomicsDBImporter {

  public static void main(final String[] args)
    throws IOException, GenomicsDBException, ParseException
  {
    if (args.length < 3) {
      System.out.println("Usage: ExampleGenomicsDBImporter" + " chromosome:interval " +
        "genomicsdbworkspace variantfile");
      System.exit(-1);
    }

    String chromosomeInterval = args[0];
    String[] temp0 = chromosomeInterval.split(":");
    String chromosomeName = temp0[0];
    String[] interval = temp0[1].split("-");
    String workspace = args[1];
    List<String> files = new ArrayList<>();
    for (int i = 2; i < args.length; ++i) {
      files.add(args[i]);
    }

    Map<String, FeatureReader<VariantContext>> map = new HashMap<>();
    List<VCFHeader> headers = new ArrayList<>();

    for (String file : files) {
      AbstractFeatureReader<VariantContext, LineIterator> reader =
        AbstractFeatureReader.getFeatureReader(file, new VCFCodec(), false);
      String sampleName = ((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0);
      headers.add((VCFHeader) reader.getHeader());
      map.put(sampleName, reader);
    }

    Set<VCFHeaderLine> mergedHeader = VCFUtils.smartMergeHeaders(headers, true);

    GenomicsDBImporter importer = new GenomicsDBImporter(
      map, mergedHeader,
      new ChromosomeInterval(chromosomeName, Integer.parseInt(interval[0]), Integer.parseInt(interval[1])),
      workspace, "", 1000L, 1048576L);
    boolean isdone = importer.importBatch();
    assert (isdone);
  }
}
