package com.intel.genomicsdb.importer.extensions;

import com.googlecode.protobuf.format.JsonFormat;
import com.intel.genomicsdb.model.GenomicsDBCallsetsMapProto;
import com.intel.genomicsdb.model.GenomicsDBImportConfiguration;
import com.intel.genomicsdb.model.GenomicsDBVidMapProto;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;

import static com.googlecode.protobuf.format.JsonFormat.printToString;

public interface JsonFileExtensions {
    /**
     * Create a JSON file from the import configuration
     *
     * @param importConfiguration The configuration object
     * @param filename            File to dump the loader JSON to
     * @return New file (with the specified name) written to local storage
     * @throws FileNotFoundException PrintWriter throws this exception when file not found
     */
    default File dumpTemporaryLoaderJSONFile(
            final GenomicsDBImportConfiguration.ImportConfiguration importConfiguration,
            final String filename) throws IOException {
        String loaderJSONString = JsonFormat.printToString(importConfiguration);

        File tempLoaderJSONFile = (filename.isEmpty()) ?
                File.createTempFile("loader_", ".json") :
                new File(filename);

        if (filename.isEmpty())
            tempLoaderJSONFile.deleteOnExit();

        PrintWriter out = new PrintWriter(tempLoaderJSONFile);
        out.println(loaderJSONString);
        out.close();
        return tempLoaderJSONFile;
    }

    /**
     * Writes a JSON file from a vidmap protobuf object. This
     * method is not called implicitly by GenomicsDB constructor.
     * It needs to be called explicitly if the user wants these objects
     * to be written. Called explicitly from GATK-4 GenomicsDBImport tool
     *
     * @param outputCallsetMapJSONFilePath Full path of file to be written
     * @param callsetMappingPB             Protobuf callset map object
     * @throws FileNotFoundException PrintWriter throws this exception when file not found
     */
    default void writeCallsetMapJSONFile(final String outputCallsetMapJSONFilePath,
                                         final GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB)
            throws FileNotFoundException {
        String callsetMapJSONString = printToString(callsetMappingPB);
        File callsetMapJSONFile = new File(outputCallsetMapJSONFilePath);
        PrintWriter out = new PrintWriter(callsetMapJSONFile);
        out.println(callsetMapJSONString);
        out.close();
    }

    /**
     * Writes a JSON file from a set of VCF header lines. This
     * method is not called implicitly by GenomicsDB constructor.
     * It needs to be called explicitly if the user wants these objects
     * to be written. Called explicitly from GATK-4 GenomicsDBImport tool
     *
     * @param outputVidMapJSONFilePath Full path of file to be written
     * @param vidMappingPB             ProtoBuf data structure for vid mapping
     * @throws FileNotFoundException PrintWriter throws this exception when file not found
     */
    default void writeVidMapJSONFile(final String outputVidMapJSONFilePath,
                                     final GenomicsDBVidMapProto.VidMappingPB vidMappingPB)
            throws FileNotFoundException {
        String vidMapJSONString = printToString(vidMappingPB);
        File vidMapJSONFile = new File(outputVidMapJSONFilePath);
        PrintWriter out = new PrintWriter(vidMapJSONFile);
        out.println(vidMapJSONString);
        out.close();
    }

    /**
     * Writes a VCF Header file from a set of VCF header lines. This
     * method is not called implicitly by GenomicsDB constructor.
     * It needs to be called explicitly if the user wants these objects
     * to be written. Called explicitly from GATK-4 GenomicsDBImport tool
     *
     * @param outputVcfHeaderFilePath Full path of file to be written
     * @param headerLines             Set of header lines
     * @throws FileNotFoundException PrintWriter throws this exception when file not found
     */
    default void writeVcfHeaderFile(final String outputVcfHeaderFilePath, final Set<VCFHeaderLine> headerLines)
            throws FileNotFoundException {
        File vcfHeaderFile = new File(outputVcfHeaderFilePath);
        final VCFHeader vcfHeader = new VCFHeader(headerLines);
        VariantContextWriter vcfWriter = new VariantContextWriterBuilder()
                .clearOptions()
                .setOutputFile(vcfHeaderFile)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .build();
        vcfWriter.writeHeader(vcfHeader);
        vcfWriter.close();
    }
}
