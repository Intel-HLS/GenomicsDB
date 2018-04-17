package com.intel.genomicsdb.reader;

import java.util.Comparator;

import static com.intel.genomicsdb.Constants.CHROMOSOME_FOLDER_DELIMITER_SYMBOL_REGEX;

public class ChrArrayFolderComparator implements Comparator<String> {
    @Override
    public int compare(String o1, String o2) {
        int chromCompare = extractChromsomeName(o1).compareTo(extractChromsomeName(o2));
        if (chromCompare == 0) {
            return extractIntervalStart(o1) - extractIntervalStart(o2);
        } else {
            return chromCompare;
        }
    }

    private String extractChromsomeName(String s) {
        String[] values = s.split(CHROMOSOME_FOLDER_DELIMITER_SYMBOL_REGEX);
        return values[0];
    }

    private int extractIntervalStart(String s) {
        String[] values = s.split(CHROMOSOME_FOLDER_DELIMITER_SYMBOL_REGEX);
        return Integer.parseInt(values[1]);
    }
}
