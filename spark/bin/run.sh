#!/bin/bash

if [ -z "$SPARK_HOME" ]; then
  echo "SPARK_HOME not set"
  exit
fi

CLASS=com.intel.genomicsdb.GenomicsDBSparkFactory
GENOMICSDB_JAR=genomicsdb-sparkapi-0.1-SNAPSHOT.jar
GENOMICSDB_SPARK_HOME="$(cd `dirname $0`; pwd)"
HTSJDK_JAR=$HOME/Downloads/htsjdk-2.4.1.jar

$SPARK_HOME/bin/spark-submit \
  --deploy-mode client \
  --class $CLASS \
  --driver-class-path $GENOMICSDB_SPARK_HOME/../target/$GENOMICSDB_JAR:$HTSJDK_JAR \
  --driver-library-path $GENOMICSDB_SPARK_HOME/../lib:/opt/gcc-5.1.0/lib64:$LD_LIBRARY_PATH \
  --conf spark.executor.extraLibraryPath=$GENOMICSDB_SPARK_HOME/lib:/opt/gcc-5.1.0/lib64:$LD_LIBRARY_PATH \
  --conf spark.executor.extraClassPath=$GENOMICSDB_SPARK_HOME/../target/$GENOMICSDB_JAR:$HTSJDK_JAR \
  --conf spark.tiledb.workspace.dir=$TILEDB_WORKSPACE \
  target/$GENOMICSDB_JAR \
  "$@"
