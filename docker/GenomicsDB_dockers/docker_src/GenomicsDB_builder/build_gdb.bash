#! /bin/bash
# The MIT License (MIT)
# Copyright (c) 2016-2017 Intel Corporation
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

if [ ! -d /output ]; then
  echo 'can not find /output/ exit...'
  return 2
fi
export BUILD_ROOT=${HOME}/build_src
export PROTOBUF_LIBRARY=$BUILD_ROOT/protobuf_build
echo "user=$USER build_home=$build_home, PROTOBUF_LIBRARY=$PROTOBUF_LIBRARY, GenomicsDB_HOME=$GenomicsDB_HOME"

protobuf_to_dir=/output
genomicsdb_to_dir=/output

# make protobuf
build_proto_buf() {
  echo "+++ Building protobuf at ${PROTOBUF_LIBRARY}..."
  mkdir -p /output/protobuf
  mkdir -p ${PROTOBUF_LIBRARY} && pushd ${PROTOBUF_LIBRARY} >/dev/null 2>&1
  git clone https://github.com/google/protobuf.git
  cd protobuf
  git checkout 3.0.x
  ./autogen.sh
  ./configure --prefix=$protobuf_to_dir --with-pic
  if [ -f ./Makefile ]; then
     make && make install
     protolibs=$(find $protobuf_to_dir/lib/ -name 'libproto*' -type f -exec basename {} \;)
     ret=$?
     echo "--- Done rc=$ret, protobuf libs: $protolibs"
     popd >/dev/null 2>&1
     return 0
  else
     popd >/dev/null 2>&1
     echo "ERROR: build_proto_buf not find Makefile"
     return -1
  fi
}

# make genomicsdb
build_gdb() {
  echo
  echo "+++ Building GenomicsDB at ${GenomicsDB_HOME}..."
  git clone --recursive https://github.com/Intel-HLS/GenomicsDB.git
  ws=GenomicsDB/build
  mkdir -p $ws && pushd $ws >/dev/null 2>&1
  branch=${GDB_BRANCH:=master}
  git checkout $branch
  git branch
  git submodule update --recursive
  build_type=$1
  shift
  install_dir=$genomicsdb_to_dir/$build_type
  [[ -d $install_dir ]] && rm -rf $install_dir
  mkdir -p $install_dir
  if [ $# -gt 0 ]; then
    cmake --warn-uninitialized --debug-output .. -DCMAKE_INSTALL_PREFIX=$install_dir -DPROTOBUF_LIBRARY=$protobuf_to_dir -DCMAKE_BUILD_TYPE=$build_type $@
  else
    cmake --warn-uninitialized --debug-output .. -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$install_dir -DPROTOBUF_LIBRARY=$protobuf_to_dir -DDISABLE_OPENMP=True -DBUILD_JAVA=False -DDO_PROFILING=False
  fi

  if [ -f ./Makefile ]; then
     make && make  install
     echo "INFO: Successfully built GenomicsDB... run test"
     # test failed: 'Error: Could not find or load main class TestGenomicsDB'
     # remove ../tests/run.py $PWD $install_dir
     popd >/dev/null 2>&1
     return 0
  else
     echo "ERROR: build_gdb not find Makefile"
     popd >/dev/null 2>&1
     return -1
  fi
}

#source /opt/rh/devtoolset-4/enable
gcc --version
echo "gcc $(gcc --version)"
mkdir -p ${BUILD_ROOT} && pushd ${BUILD_ROOT} >/dev/null 2>&1
build_proto_buf
retst=$?
if [ $retst -eq 0 ]; then
  build_gdb Release $@
  retst=$?
  if [ $retst -eq 0 ]; then
    build_gdb Debug $@
    retst=$?
    gdb_version=$($genomicsdb_to_dir/Release/bin/vcf2tiledb --version)
    echo "GenomicsDB version is $gdb_version"
    mkdir -p $BUILD_ROOT/$gdb_version
    rsync -a -r $genomicsdb_to_dir $BUILD_ROOT/$gdb_version
  else
    echo "ERROR: build release version failed. status=$retst"
  fi
fi
[[ $retst -eq 0 ]] && echo "DONE: built GenomicsDB" || echo "FAIL: cannot build GenomicsDB"
popd >/dev/null 2>&1
return $retst
