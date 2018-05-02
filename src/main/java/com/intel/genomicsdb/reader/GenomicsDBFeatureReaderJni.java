package com.intel.genomicsdb.reader;

class GenomicsDBFeatureReaderJni {

    /**
     * Checks if GenomicsDB array exists.
     * @param workspace
     * @param arrayName
     * @return <code>true</code> if workspace with arrayName exists else return <code>false</code>
     */
    public static native boolean jniIsTileDBArray(final String workspace, final String arrayName);
}
