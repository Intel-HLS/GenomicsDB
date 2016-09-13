package com.intel.genomicsdb;

import org.apache.hadoop.conf.Configuration;
import org.apache.log4j.Logger;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class GenomicsDBConf extends Configuration implements Serializable {

  static final String LOADERJSON = "genomicsdb.input.loaderjsonfile";
  static final String QUERYJSON = "genomicsdb.input.queryjsonfile";
  static final String MPIHOSTFILE = "genomicsdb.input.mpi.hostfile";

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

  public GenomicsDBConf setHostFile(String path) throws FileNotFoundException {
    set(MPIHOSTFILE, path);

    Logger logger = Logger.getLogger(GenomicsDBConf.class);
    Scanner scanner = new Scanner(new FileInputStream(path));
    while (scanner.hasNextLine()) {
      String host = scanner.nextLine();
      hosts.add(host);
      logger.error("host file content: " + host);
    }
    return this;
  }

  List<String> getHosts() {
    return hosts;
  }

  public Configuration value() {
    return this;
  }
}

