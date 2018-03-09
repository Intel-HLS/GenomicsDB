#!/bin/bash

if [ -z "$SPARK_HOME" ]; then
  echo "SPARK_HOME not set"
  exit
fi

SPARK_MASTER_URL="$1"
shift

CLASS=com.intel.genomicsdb.spark.GenomicsDBJavaSparkFactory
GENOMICSDB_JAR=genomicsdb-0.4.0-jar-with-dependencies.jar
GENOMICSDB_BIN_DIR="$(cd `dirname $0`; pwd)"
GENOMICSDB_JAR_PATH=$GENOMICSDB_BIN_DIR/$GENOMICSDB_JAR

CLASSPATH=$CLASSPATH:$GENOMICSDB_JAR_PATH

$SPARK_HOME/bin/spark-submit \
  --class $CLASS \
  --master $SPARK_MASTER_URL \
  --deploy-mode client \
  "$GENOMICSDB_JAR_PATH" \
  "$@"
