package com.intel.genomicsdb.reader;


import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class ChrArrayFolderComparatorSpec {
    @Test(testName = "should compare first by chromosome and then interval start")
    public void shouldCompareFirstByChromosomeAndThenIntervalStart() {
        //Given
        List<String> arrayList = new ArrayList<>();
        arrayList.add("1#100#200");
        arrayList.add("1#201#300");
        arrayList.add("2#100#200");
        arrayList.add("2#201#300");
        arrayList.add("1#1000#2000");

        //When
        arrayList.sort(new ChrArrayFolderComparator());

        //Then
        Assert.assertEquals(arrayList.get(2), "1#1000#2000");
    }
}
