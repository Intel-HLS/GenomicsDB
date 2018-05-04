package com.intel.genomicsdb.reader;


import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

import static com.intel.genomicsdb.Constants.CHROMOSOME_INTERVAL_FOLDER;

public class ChrArrayFolderComparatorSpec {
    @Test(testName = "should compare first by chromosome and then interval start")
    public void shouldCompareFirstByChromosomeAndThenIntervalStart() {
        //Given
        List<String> arrayList = new ArrayList<>();
        arrayList.add(String.format(CHROMOSOME_INTERVAL_FOLDER, "1", 100, 200));
        arrayList.add(String.format(CHROMOSOME_INTERVAL_FOLDER, "1", 201, 300));
        arrayList.add(String.format(CHROMOSOME_INTERVAL_FOLDER, "2", 100, 200));
        arrayList.add(String.format(CHROMOSOME_INTERVAL_FOLDER, "2", 201, 300));
        arrayList.add(String.format(CHROMOSOME_INTERVAL_FOLDER, "1", 1000, 2000));

        //When
        arrayList.sort(new ChrArrayFolderComparator());

        //Then
        Assert.assertEquals(arrayList.get(2), String.format(CHROMOSOME_INTERVAL_FOLDER, "1", 1000, 2000));
    }
}
