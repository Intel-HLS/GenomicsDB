package com.intel.genomicsdb.importer;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

public final class Constants {

    private Constants() {
    }

    public static final int DEFAULT_ZERO_BATCH_SIZE = 0;
    //Allele specific annotation fields
    public static final HashSet<String> R_LENGTH_HISTOGRAM_FIELDS_FLOAT_BINS = new HashSet<>(Arrays.asList(
            "AS_RAW_BaseQRankSum",
            "AS_RAW_MQRankSum",
            "AS_RAW_ReadPosRankSum"
    ));
    public static final HashSet<String> R_LENGTH_TWO_DIM_FLOAT_VECTOR_FIELDS = new HashSet<>(Collections.singletonList(
            "AS_RAW_MQ"
    ));
    public static final HashSet<String> R_LENGTH_TWO_DIM_INT_VECTOR_FIELDS = new HashSet<>(Collections.singletonList(
            "AS_SB_TABLE"
    ));
    public static final long DEFAULT_BUFFER_CAPACITY = 20480; //20KB
}