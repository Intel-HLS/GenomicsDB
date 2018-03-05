/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

package com.intel.genomicsdb;
import java.io.PrintStream;
import java.lang.management.ManagementFactory;

public class GenomicsDBTimer
{
    private long mBeginWallClockTime = 0;
    private long mBeginCpuTime = 0;

    //Interval time
    private long mLastIntervalWallClockTime = 0;
    private long mLastIntervalCpuTime = 0;

    //Cumulative time
    private double mCumulativeWallClockTime = 0;
    private double mCumulativeCpuTime = 0;

    public GenomicsDBTimer()
    {
        start();
    }

    public void start()
    {
        mBeginWallClockTime = System.nanoTime();
        mBeginCpuTime = ManagementFactory.getThreadMXBean().getCurrentThreadCpuTime();
    }

    public void stop()
    {
        long endWallClockTime = System.nanoTime();
        long endCpuTime = ManagementFactory.getThreadMXBean().getCurrentThreadCpuTime();
        mLastIntervalWallClockTime = (endWallClockTime - mBeginWallClockTime);
        mLastIntervalCpuTime = (endCpuTime - mBeginCpuTime);
        mCumulativeWallClockTime += (((double)mLastIntervalWallClockTime)/1e9);
        mCumulativeCpuTime += (((double)mLastIntervalCpuTime)/1e9);
    }

    public void print(final String prefix, PrintStream fptr)
    {
        fptr.print("GENOMICSDB_TIMER,");
        if(!prefix.isEmpty())
            fptr.print(prefix+',');
        fptr.println("Wall-clock time(s)," + mCumulativeWallClockTime + ",Cpu time(s),"
                + mCumulativeCpuTime);
    }
}
