#
# FindHDFS.cmake
#
#
# The MIT License
#
# Copyright (c) 2017 Omics Data Automation, Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

# Get java and hdfs include directories and libraries for linking in hdfs support
# This module sets the following result variables
#  For include
#   JNI_INCLUDE_DIRS
#   HDFS_INCLUDE_PATH
#  For link
#   JAVA_JVM_LIBRARY
#   HDFS_LIBRARY

# Using system provided FindJNI.cmake
find_package(JNI)
if (NOT JNI_FOUND)
  message(FATAL_ERROR "JAVA installation not found")
  return()
endif()
message(STATUS "Found JVM: " ${JAVA_JVM_LIBRARY}) 

set(HADOOP_HOME "$ENV{HADOOP_HOME}")
if(NOT HADOOP_HOME)
  message(FATAL_ERROR "Undefined HADOOP_HOME environment variable")
  return()
endif()
if(NOT EXISTS ${HADOOP_HOME})
  message(FATAL_ERROR "HADOOP_HOME=" ${HADOOP_HOME} " NOT found")
  return()
endif()
if (NOT IS_DIRECTORY ${HADOOP_HOME})
  message(FATAL_ERROR "HADOOP_HOME=" ${HADOOP_HOME} " NOT a directory")
  return()
endif()

message(STATUS "Found env HADOOP_HOME: " ${HADOOP_HOME})

if (IS_DIRECTORY ${HADOOP_HOME}/include)
  set(HDFS_INCLUDE_PATH ${JNI_INCLUDE_DIRS} ${HADOOP_HOME}/include)
else()
  message(FATAL_ERROR "HADOOP_HOME=" ${HADOOP_HOME}
	" does NOT seem to contain an include folder")
  return()
endif()

find_library(HDFS_LIBRARY hdfs PATHS ${HADOOP_HOME}/lib/native)
if(HDFS_LIBRARY)
  message(STATUS "Found hdfs: " ${HDFS_LIBRARY})
else()
  message(FATAL_ERROR "hdfs library could not be found")
  return()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDFS "Could not find libhdfs libraries ${DEFAULT_MSG}" HDFS_INCLUDE_PATH HDFS_LIBRARY)

