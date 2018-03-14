package com.intel.genomicsdb.importer.extensions;

import com.intel.genomicsdb.model.GenomicsDBVidMapProto;
import htsjdk.variant.vcf.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import static com.intel.genomicsdb.importer.Constants.R_LENGTH_HISTOGRAM_FIELDS_FLOAT_BINS;
import static com.intel.genomicsdb.importer.Constants.R_LENGTH_TWO_DIM_FLOAT_VECTOR_FIELDS;
import static com.intel.genomicsdb.importer.Constants.R_LENGTH_TWO_DIM_INT_VECTOR_FIELDS;

public interface VidMapExtensions {
    /**
     * Generate the ProtoBuf data structure for vid mapping
     * Also remember, contigs are 1-based which means
     * TileDB column offsets should start from 1
     *
     * @param mergedHeader Header from all input GVCFs are merged to
     *                     create one vid map for all
     * @return a vid map containing all field names, lengths and types
     * from the merged GVCF header
     */
    default GenomicsDBVidMapProto.VidMappingPB generateVidMapFromMergedHeader(final Set<VCFHeaderLine> mergedHeader) {
        List<GenomicsDBVidMapProto.InfoField> infoFields = new ArrayList<>();
        List<GenomicsDBVidMapProto.Chromosome> contigs = new ArrayList<>();

        //ID field
        GenomicsDBVidMapProto.InfoField.Builder IDFieldBuilder =
                GenomicsDBVidMapProto.InfoField.newBuilder();
        GenomicsDBVidMapProto.FieldLengthDescriptorComponentPB.Builder lengthDescriptorComponentBuilder =
                GenomicsDBVidMapProto.FieldLengthDescriptorComponentPB.newBuilder();
        lengthDescriptorComponentBuilder.setVariableLengthDescriptor("var");
        IDFieldBuilder.setName("ID").addType("char").addLength(lengthDescriptorComponentBuilder.build());
        infoFields.add(IDFieldBuilder.build());

        int dpIndex = -1;
        long columnOffset = 0L;

        for (VCFHeaderLine headerLine : mergedHeader) {
            GenomicsDBVidMapProto.InfoField.Builder infoBuilder =
                    GenomicsDBVidMapProto.InfoField.newBuilder();

            if (headerLine instanceof VCFFormatHeaderLine) {
                VCFFormatHeaderLine formatHeaderLine = (VCFFormatHeaderLine) headerLine;
                boolean isGT = formatHeaderLine.getID().equals(VCFConstants.GENOTYPE_KEY);
                String genomicsDBType = isGT ? "int" : formatHeaderLine.getType().toString();
                String genomicsDBLength = isGT ? "PP" : (formatHeaderLine.getType() == VCFHeaderLineType.String)
                        ? "VAR" : getVcfHeaderLength(formatHeaderLine);
                lengthDescriptorComponentBuilder.setVariableLengthDescriptor(genomicsDBLength);

                infoBuilder
                        .setName(formatHeaderLine.getID())
                        .addType(genomicsDBType)
                        .addLength(lengthDescriptorComponentBuilder.build());

                if (formatHeaderLine.getID().equals("DP") && dpIndex != -1) {
                    GenomicsDBVidMapProto.InfoField prevDPField = removeInfoFields(infoFields, dpIndex);

                    infoBuilder
                            .addVcfFieldClass(prevDPField.getVcfFieldClass(0))
                            .addVcfFieldClass("FORMAT");
                } else {
                    infoBuilder
                            .addVcfFieldClass("FORMAT");
                }

                GenomicsDBVidMapProto.InfoField formatField = infoBuilder.build();
                infoFields.add(formatField);
                if (formatHeaderLine.getID().equals("DP")) {
                    dpIndex = infoFields.indexOf(formatField);
                }
            } else if (headerLine instanceof VCFInfoHeaderLine) {
                VCFInfoHeaderLine infoHeaderLine = (VCFInfoHeaderLine) headerLine;
                final String infoFieldName = infoHeaderLine.getID();
                infoBuilder
                        .setName(infoFieldName);
                //allele specific annotations
                if (R_LENGTH_HISTOGRAM_FIELDS_FLOAT_BINS.contains(infoFieldName)
                        || R_LENGTH_TWO_DIM_FLOAT_VECTOR_FIELDS.contains(infoFieldName)
                        || R_LENGTH_TWO_DIM_INT_VECTOR_FIELDS.contains(infoFieldName)
                        ) {
                    lengthDescriptorComponentBuilder.setVariableLengthDescriptor("R");
                    infoBuilder.addLength(lengthDescriptorComponentBuilder.build());
                    lengthDescriptorComponentBuilder.setVariableLengthDescriptor("var"); //ignored - can set anything here
                    infoBuilder.addLength(lengthDescriptorComponentBuilder.build());
                    infoBuilder.addVcfDelimiter("|");
                    infoBuilder.addVcfDelimiter(",");
                    if (R_LENGTH_HISTOGRAM_FIELDS_FLOAT_BINS.contains(infoFieldName)) {
                        //Each element of the vector is a tuple <float, int>
                        infoBuilder.addType("float");
                        infoBuilder.addType("int");
                        infoBuilder.setVCFFieldCombineOperation("histogram_sum");
                    } else {
                        infoBuilder.setVCFFieldCombineOperation("element_wise_sum");
                        if (R_LENGTH_TWO_DIM_FLOAT_VECTOR_FIELDS.contains(infoFieldName)) {
                            infoBuilder.addType("float");
                        } else if (R_LENGTH_TWO_DIM_INT_VECTOR_FIELDS.contains(infoFieldName)) {
                            infoBuilder.addType("int");
                        }
                    }
                } else {
                    lengthDescriptorComponentBuilder.setVariableLengthDescriptor(
                            infoHeaderLine.getType() == VCFHeaderLineType.String ? "var" : getVcfHeaderLength(infoHeaderLine));

                    infoBuilder.addType(infoHeaderLine.getType().toString())
                            .addLength(lengthDescriptorComponentBuilder.build());
                }

                if (infoHeaderLine.getID().equals("DP") && dpIndex != -1) {
                    GenomicsDBVidMapProto.InfoField prevDPield = removeInfoFields(infoFields, dpIndex);

                    infoBuilder
                            .addVcfFieldClass(prevDPield.getVcfFieldClass(0))
                            .addVcfFieldClass("INFO");
                } else {
                    infoBuilder
                            .addVcfFieldClass("INFO");
                }

                GenomicsDBVidMapProto.InfoField infoField = infoBuilder.build();
                infoFields.add(infoField);
                if (infoHeaderLine.getID().equals("DP")) {
                    dpIndex = infoFields.indexOf(infoField);
                }
            } else if (headerLine instanceof VCFFilterHeaderLine) {
                VCFFilterHeaderLine filterHeaderLine = (VCFFilterHeaderLine) headerLine;

                infoBuilder.setName(filterHeaderLine.getID());

                if (!filterHeaderLine.getValue().isEmpty()) {
                    infoBuilder.addType(filterHeaderLine.getValue());
                } else {
                    infoBuilder.addType("int");
                }
                infoBuilder.addVcfFieldClass("FILTER");
                GenomicsDBVidMapProto.InfoField filterField = infoBuilder.build();

                infoFields.add(filterField);
            } else if (headerLine instanceof VCFContigHeaderLine) {
                VCFContigHeaderLine contigHeaderLine = (VCFContigHeaderLine) headerLine;

                long length = contigHeaderLine.getSAMSequenceRecord().getSequenceLength();
                GenomicsDBVidMapProto.Chromosome.Builder contigBuilder =
                        GenomicsDBVidMapProto.Chromosome.newBuilder();
                contigBuilder
                        .setName(contigHeaderLine.getID())
                        .setLength(length)
                        .setTiledbColumnOffset(columnOffset);

                columnOffset += length;

                GenomicsDBVidMapProto.Chromosome chromosome = contigBuilder.build();

                contigs.add(chromosome);
            }
        }

        GenomicsDBVidMapProto.VidMappingPB.Builder vidMapBuilder =
                GenomicsDBVidMapProto.VidMappingPB.newBuilder();

        return vidMapBuilder
                .addAllFields(infoFields)
                .addAllContigs(contigs)
                .build();
    }

    /**
     * Maps the "Number" from INFO or FORMAT fields in VCF header
     * to GenomicsDB lengths. If unbounded, the field becomes a
     * variable length attribute (discussed in more detail in TileDB
     * tutorial tiledb.org). If "A", "R" or "G", the field also
     * becomes a variable length attribute, but the order and length
     * depends on order/number of alleles or genotypes. Integers
     * remain integers
     *
     * @param headerLine Info or Format header line from VCF
     * @return VAR, A, R, G, or integer values from VCF
     */
    default String getVcfHeaderLength(final VCFHeaderLine headerLine) {
        VCFHeaderLineCount type;

        if (headerLine instanceof VCFFormatHeaderLine) {
            type = ((VCFFormatHeaderLine) headerLine).getCountType();
        } else {
            type = ((VCFInfoHeaderLine) headerLine).getCountType();
        }

        String length = "";
        int count;
        boolean isFlagType;
        switch (type) {
            case UNBOUNDED:
                length = "VAR";
                break;
            case A:
                length = "A";
                break;
            case R:
                length = "R";
                break;
            case G:
                length = "G";
                break;
            case INTEGER: {
                if (headerLine instanceof VCFFormatHeaderLine) {
                    VCFFormatHeaderLine formatHeaderLine = (VCFFormatHeaderLine) headerLine;
                    count = formatHeaderLine.getCount();
                    isFlagType = formatHeaderLine.getType().equals(VCFHeaderLineType.Flag);
                } else {
                    VCFInfoHeaderLine infoHeaderLine = (VCFInfoHeaderLine) headerLine;
                    count = infoHeaderLine.getCount();
                    isFlagType = infoHeaderLine.getType().equals(VCFHeaderLineType.Flag);
                }
                //Weird Flag fields - Number=0 in the VCF header :(
                if (count == 0 && isFlagType)
                    length = "1";
                else
                    length = String.valueOf(count);
                break;
            }
        }
        return length;
    }

    default GenomicsDBVidMapProto.InfoField removeInfoFields(List<GenomicsDBVidMapProto.InfoField> infoFields,
                                                             final int dpIndex) {
        GenomicsDBVidMapProto.InfoField dpFormatField = infoFields.get(dpIndex);
        infoFields.remove(dpIndex);
        return dpFormatField;
    }
}
