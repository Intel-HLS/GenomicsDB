package com.intel.genomicsdb.reader;

class GenomicsDBFeatureReaderJni {

    /**
     * Checks if GenomicsDB array exists.
     * @param workspace workspace
     * @param arrayName array name
     * @return <code>true</code> if workspace with arrayName exists else return <code>false</code>
     */
    public static native boolean jniIsTileDBArray(final String workspace, final String arrayName);

    /**
     * List Arrays in given workspace.
     * @param workspace workspace
     * @return list of arrays in given workspace.
     */
    public static native String[] jniListTileDBArrays(final String workspace);
}
