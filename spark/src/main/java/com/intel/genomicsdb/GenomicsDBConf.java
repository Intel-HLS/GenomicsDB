package com.intel.genomicsdb;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;

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
    setLoaderJsonFile(new Path(configuration.get(GenomicsDBConf.LOADERJSON)));
    setQueryJsonFile(new Path(configuration.get(GenomicsDBConf.QUERYJSON)));
    setHostFile(new Path(configuration.get(GenomicsDBConf.MPIHOSTFILE)));
  }

  // <String> left for backward compatibility to Java 7
  private ArrayList<String> hosts = new ArrayList<>();

//  public GenomicsDBConf setMinInputSplitSize(Job job, Long size) {
//    job.getConfiguration().setLong(SPLIT_MINSIZE, size);
//    return this;
//  }

//  public Long getMinSplitSize(JobContext job) {
//    return job.getConfiguration().getLong(SPLIT_MINSIZE, 1L);
//  }

  public GenomicsDBConf setLoaderJsonFile(Path path) {
    set(LOADERJSON, path.toString());
    return this;
  }

//  public Path getLoaderJsonFile() {
//    return new Path(get(LOADERJSON));
//  }

  public GenomicsDBConf setQueryJsonFile(Path path) {
    set(QUERYJSON, path.toString());
    return this;
  }

//  public Path getQueryJsonFile() {
//    return new Path(get(QUERYJSON));
//  }

  public GenomicsDBConf setHostFile(Path path) throws FileNotFoundException {
    set(MPIHOSTFILE, path.toString());

    Scanner scanner = new Scanner(new FileInputStream(path.toString()));
    while (scanner.hasNextLine()) {
      hosts.add(scanner.nextLine());
    }
    return this;
  }

//  public Path getHostFile(){
//    return new Path(get(MPIHOSTFILE));
//  }

  List<String> getHosts() {
    return hosts;
  }

  public Configuration value() {
    return this;
  }
}

