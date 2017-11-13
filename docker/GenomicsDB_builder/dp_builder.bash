#! /bin/bash

# License and Copyright
# GenomicsDB Builder Docker was created while Ming Rutar worked at Intel. The code is not
# part of assignment. GenomicsDB is an open source project, see https://github.com/Intel-HLS/GenomicsDB
# The project is stored under Ming's Intel email account and subject to Intel copyright.
# Included GenomicsDB Intel copyright below.
# 
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

# $MY_DOCKER_ROOT directory of the project
DOCKER_ROOT=${MY_DOCKER_ROOT}

echo "START building GenomicsDB builder docker...."
pushd $DOCKER_ROOT/GDB-builder_centos >/dev/null 2>&1
$(docker volume inspect -f "{{json .Mountpoint}}" $GENOME_VOLUME) >/dev/null 2>&1
if [ $? -ne 0 ]; then
    docker volume create --name $GENOME_VOLUME
    echo "created volume $GENOME_VOLUME, rc=$?"
fi
cp -f $DOCKER_ROOT/build_gdb.bash ./usr/bin/build_genomicsdb
chmod +x ./usr/bin/build_genomicsdb
[[ -d ./usr/share/cont-layer/common/env ]] || mkdir -p ./usr/share/cont-layer/common/env/
cp -f $DOCKER_ROOT/enabledevtoolset-4.sh ./usr/share/cont-layer/common/env/enabledevtoolset-4.sh
chmod +x ./usr/share/cont-layer/common/env/enabledevtoolset-4.sh

# if you do not use proxy, remove the build args 
docker build -t genomicsdb_builder --no-cache --build-arg http_proxy=$http_proxy --build-arg https_proxy=$http_proxy . 
echo "DONE building docker....rc=$?"
# run
popd >/dev/null 2>&1
