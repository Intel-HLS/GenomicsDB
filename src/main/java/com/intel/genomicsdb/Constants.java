package com.intel.genomicsdb;

public final class Constants {

    private Constants() {
    }

    public static final String CHROMOSOME_FOLDER_DELIMITER_SYMBOL = "$";
    public static final String CHROMOSOME_FOLDER_DELIMITER_SYMBOL_REGEX = "\\$";
    public static final String CHROMOSOME_INTERVAL_FOLDER = String.format("%%s%s%%d%s%%d",
            CHROMOSOME_FOLDER_DELIMITER_SYMBOL, CHROMOSOME_FOLDER_DELIMITER_SYMBOL);
}
