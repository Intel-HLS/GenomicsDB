/*
 * The MIT License (MIT)
 * Copyright (c) 2016 Intel Corporation
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

import org.apache.hadoop.conf.Configuration;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class GenomicsDBConf extends Configuration implements Serializable {

  public static final String LOADERJSON = "genomicsdb.input.loaderjsonfile";
  public static final String QUERYJSON = "genomicsdb.input.queryjsonfile";
  public static final String MPIHOSTFILE = "genomicsdb.input.mpi.hostfile";

  public GenomicsDBConf(Configuration configuration) throws FileNotFoundException {
    super(configuration);
  }

  // <String> left for backward compatibility to Java 7
  private ArrayList<String> hosts = new ArrayList<>();

  public GenomicsDBConf setLoaderJsonFile(String path) {
    set(LOADERJSON, path);
    return this;
  }

  public GenomicsDBConf setQueryJsonFile(String path) {
    set(QUERYJSON, path);
    return this;
  }

  /**
   * Host file contains the hosts where GenomicsDB instances reside.
   * This file can be a replica of the slaves file. For now, we have
   * kept it separate, can be merged later
   *
   * @param path  Full path of the host file
   * @return  GenomicsDBConf object
   * @throws FileNotFoundException  If file not found, throw exception
   */
  public GenomicsDBConf setHostFile(String path) throws FileNotFoundException {
    set(MPIHOSTFILE, path);

    Scanner scanner = new Scanner(new FileInputStream(path));
    while (scanner.hasNextLine()) {
      String host = scanner.nextLine();
      hosts.add(host);
    }
    return this;
  }

  List<String> getHosts() {
    return hosts;
  }
}

