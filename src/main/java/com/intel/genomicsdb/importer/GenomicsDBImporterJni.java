package com.intel.genomicsdb.importer;

class GenomicsDBImporterJni {
    /**
     * Obtain the chromosome intervals for the column partition specified in the loader JSON file
     * identified by the rank. The information is returned as a string in JSON format
     * {
     * "contigs": [
     * { "chr1": [ 100, 200] },
     * { "chr2": [ 500, 600] }
     * ]
     * }
     *
     * @param loaderJSONFile path to loader JSON file
     * @param rank           rank/partition index
     * @return chromosome intervals for the queried column partition in JSON format
     */
    static native String jniGetChromosomeIntervalsForColumnPartition(final String loaderJSONFile, final int rank);

    /**
     * Create TileDB workspace
     *
     * @param workspace path to workspace directory
     * @param overwrite if set existing directory is deleted
     * @return status 0 = workspace created,
     * -1 = path was not a directory,
     * -2 = failed to create workspace,
     * 1 = existing directory, nothing changed
     */
    static native int jniCreateTileDBWorkspace(final String workspace, final boolean overwrite);

    /**
     * Write contents into file
     * @param filename path to file
     * @param contents buffer to be written out
     * @param length of buffer to be written out
     * @return status 0 = OK
     */
    static native int jniWriteToFile(final String filename, final String contents, final long length);

    /**
     * Copy source path contents to destination
     * @param source path to source file
     * @param destination path to destination file
     * @return status 0 = OK
     */
    static native int jniMoveFile(final String source, final String destination);

    /**
     * Consolidate TileDB array
     *
     * @param workspace path to workspace directory
     * @param arrayName array name
     */
    static native void jniConsolidateTileDBArray(final String workspace, final String arrayName);

    /**
     * Creates GenomicsDBImporter object when importing VCF files (no streams)
     *
     * @param loaderJSONFile Path to loader JSON file
     * @param rank           Rank of object - corresponds to the partition index in the loader
     *                       for which this object will import data
     * @param lbRowIdx       Smallest row idx which should be imported by this object
     * @param ubRowIdx       Largest row idx which should be imported by this object
     * @return status - 0 if everything was ok, -1 otherwise
     */
    native int jniGenomicsDBImporter(String loaderJSONFile, int rank, long lbRowIdx, long ubRowIdx);

    /**
     * Creates GenomicsDBImporter object when importing VCF files (no streams)
     *
     * @param loaderJSONFile Path to loader JSON file
     * @param rank           Rank of object - corresponds to the partition index in the
     *                       loader for which this object will import data
     * @param lbRowIdx       Smallest row idx which should be imported by this object
     * @param ubRowIdx       Largest row idx which should be imported by this object
     * @return "pointer"/address to GenomicsDBImporter object in memory,
     * if 0, then something went wrong
     */
    native long jniInitializeGenomicsDBImporterObject(String loaderJSONFile, int rank, long lbRowIdx, long ubRowIdx);

    /**
     * Copy the vid map protocol buffer to C++ through JNI
     *
     * @param genomicsDBImporterHandle Reference to a C++ GenomicsDBImporter object
     * @param vidMapAsByteArray        INFO, FORMAT, FILTER header lines and contig positions
     * @return Reference to a C++ GenomicsDBImporter object as long
     */
    native long jniCopyVidMap(long genomicsDBImporterHandle, byte[] vidMapAsByteArray);

    /**
     * Copy the callset map protocol buffer to C++ through JNI
     *
     * @param genomicsDBImporterHandle Reference to a C++ GenomicsDBImporter object
     * @param callsetMapAsByteArray    Callset name and row index map
     * @return Reference to a C++ GenomicsDBImporter object as long
     */
    native long jniCopyCallsetMap(long genomicsDBImporterHandle, byte[] callsetMapAsByteArray);

    /**
     * Notify importer object that a new stream is to be added
     *
     * @param genomicsDBImporterHandle "pointer" returned by jniInitializeGenomicsDBImporterObject
     * @param streamName               name of the stream
     * @param isBCF                    use BCF format to pass data to C++ layer
     * @param bufferCapacity           in bytes
     * @param buffer                   initialization buffer containing the VCF/BCF header
     * @param numValidBytesInBuffer    num valid bytes in the buffer (length of the header)
     */
    native void jniAddBufferStream(long genomicsDBImporterHandle, String streamName, boolean isBCF, long bufferCapacity,
                                   byte[] buffer, long numValidBytesInBuffer);

    /**
     * Setup loader after all the buffer streams are added
     *
     * @param genomicsDBImporterHandle "pointer" returned by jniInitializeGenomicsDBImporterObject
     * @param callsetMappingJSON       JSON formatted string containing globally consistent callset
     *                                 name to row index mapping
     * @param usingVidMappingProtoBuf  use protocol buffer based header
     * @return maximum number of buffer stream identifiers that can be returned in
     * mExhaustedBufferStreamIdentifiers later
     * (this depends on the number of partitions and the number of buffer streams)
     */
    native long jniSetupGenomicsDBLoader(long genomicsDBImporterHandle, final String callsetMappingJSON,
                                         boolean usingVidMappingProtoBuf);

    /**
     * @param handle                "pointer" returned by jniInitializeGenomicsDBImporterObject
     * @param streamIdx             stream index
     * @param partitionIdx          partition index (unused now)
     * @param buffer                buffer containing data
     * @param numValidBytesInBuffer num valid bytes in the buffer
     */
    native void jniWriteDataToBufferStream(long handle, int streamIdx, int partitionIdx, byte[] buffer,
                                           long numValidBytesInBuffer);

    /**
     * Import the next batch of data into TileDB/GenomicsDB
     *
     * @param genomicsDBImporterHandle   "pointer" returned by jniInitializeGenomicsDBImporterObject
     * @param exhaustedBufferIdentifiers contains the list of exhausted buffer stream identifiers
     *                                   - the number of
     *                                   exhausted streams is stored in the last element of the array
     * @return true if the whole import process is completed, false otherwise
     */
    native boolean jniImportBatch(long genomicsDBImporterHandle, long[] exhaustedBufferIdentifiers);
}
