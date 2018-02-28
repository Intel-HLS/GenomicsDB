/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2018 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package com.intel.genomicsdb.model;


import org.testng.Assert;
import org.testng.annotations.Test;

public class CommandLineImportConfigSpec {

    @Test(testName = "should throw an exception when there are intersections within chromosome intervals",
            expectedExceptions = IllegalArgumentException.class,
            expectedExceptionsMessageRegExp = "There are multiple intervals sharing same value. This is not allowed. " +
            "Intervals should be defined without intersections.")
    public void shouldThrowExceptionWhenThereAreIntersectionsInChromosomeIntervals() {
        //Given
        String[] args = ("-L 1:12000-13000 -L 1:17000-18000 -L 1:18000-19000 -w path/to/ws " +
                "--vidmap-output path/to/vidmap --callset-output path/to/callset tests/inputs/vcfs/t0.vcf.gz").split(" ");
        //When
        new CommandLineImportConfig("TestGenomicsDBImporterWithMergedVCFHeader", args);

        //Then
        //Exception is expected
    }

    @Test(testName = "should not throw exception when there are intersections on intervals ranges but not for the same chromosome")
    public void shouldNotThrowExceptionWhenThereAreIntersectionsButDifferentChromosomes() {
        //Given
        String[] args = ("-L 1:12000-13000 -L 2:12000-13000 -w path/to/ws " +
                "--vidmap-output path/to/vidmap --callset-output path/to/callset tests/inputs/vcfs/t0.vcf.gz").split(" ");
        //When
        CommandLineImportConfig config = new CommandLineImportConfig("TestGenomicsDBImporterWithMergedVCFHeader", args);

        //Then
        Assert.assertEquals(config.getSampleNameToVcfPath().values().toArray()[0].toString(), "tests/inputs/vcfs/t0.vcf.gz");
    }

    @Test(testName = "should not throw exception when there is only one interval")
    public void shouldNotThrowExceptionWhenThereIsOnlyOneInteval() {
        //Given
        String[] args = ("-L 1:12000-13000 -w path/to/ws --vidmap-output path/to/vidmap " +
                "--callset-output path/to/callset tests/inputs/vcfs/t0.vcf.gz").split(" ");
        //When
        CommandLineImportConfig config = new CommandLineImportConfig("TestGenomicsDBImporterWithMergedVCFHeader", args);

        //Then
        Assert.assertEquals(config.getSampleNameToVcfPath().values().toArray()[0].toString(), "tests/inputs/vcfs/t0.vcf.gz");
    }
}
