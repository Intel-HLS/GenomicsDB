#!/bin/bash

if [ -z "$SPARK_HOME" ]; then
  echo "SPARK_HOME not set"
  exit
fi

if [ -z "$TILEDB_WORKSPACE" ]; then
  echo "TILEDB_WORKSPACE not set"
  exit
fi

CLASS=com.intel.genomicsdb.GenomicsDBJavaSparkFactory
GENOMICSDB_JAR=genomicsdb-0.4.0.jar
GENOMICSDB_SPARK_HOME="$(cd `dirname $0`; pwd)"

# Add full path of the htsjdk jar file here
# e.g. HTSJDK_JAR=$HOME/Downloads/htsjdk-2.4.1.jar
HTSJDK_JAR=

# Add full path of the breeze math and core jar files here
# e.g. BREEZE_JAR=$HOME/Downloads/breeze-math_2.10-0.4.jar:$HOME/Downloads/breeze-core_2.10-0.2.jar
BREEZE_JAR=

CLASSPATH+=$HTSJDK_JAR:$BREEZE_JAR:$GENOMICSDB_SPARK_HOME/../target/$GENOMICSDB_JAR

$SPARK_HOME/bin/spark-submit \
  --deploy-mode client \
  --class $CLASS \
  --driver-class-path $CLASSPATH \
  --driver-library-path $LD_LIBRARY_PATH \
  --conf spark.executor.extraLibraryPath=$LD_LIBRARY_PATH \
  --conf spark.executor.extraClassPath=$CLASSPATH \
  --conf spark.tiledb.workspace.dir=$TILEDB_WORKSPACE \
  target/$GENOMICSDB_JAR \
  "$@"
