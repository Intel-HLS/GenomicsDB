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

#! /bin/bash
# usage: $0 [GenomicsDB_branch_name]
#   use master branch if no GenomicsDB_branch_name provided

echo "START building GenomicsDB, $#.... "
if [ $# -lt 1 ]; then
    my_branch="master"
else
    my_branch="$1"
fi
host_dest_dir="${HOME}/docker-geno_${my_branch}_$(date +'%y%m%d')"

# ${GENOME_VOLUME} name of docker volume
volume_dir=$(docker volume inspect -f "{{json .Mountpoint}}" ${GENOME_VOLUME} | sed "s/^\"\(.*\)\"$/\1/")
if sudo test -d ${volume_dir}; then
  output=${GENOME_VOLUME}
else
  output="$HOME/docker_build"
  if sudo test -d $output; then
    sudo rm -rf $output
  fi
  mkdir -p $output
fi
cntr_name="GDB-${my_branch}_$(date +'%y%m%d_%H%M')"
echo "RUN my_branch is ${my_branch}, volume_dir is ${volume_dir}, host_dest_dir is ${host_dest_dir}, output is ${output} cntr_name is ${cntr_name}"
docker run -it --name ${cntr_name} -v $output:/output/ -e http_proxy=$http_proxy -e https_proxy=$http_proxy -e HTTP_PROXY=$http_proxy -e HTTPS_PROXY=$http_proxy -e GDB_BRANCH=${my_branch} genomicsdb_builder build_genomicsdb
echo "DONE building GenomicsDB, output at $output.... start copy "

sudo rsync -av ${volume_dir}/ ${host_dest_dir}
sudo chown -R $(whoami): ${host_dest_dir}
ln -s ${host_dest_dir}/Release $(dirname ${host_dest_dir})/docker-geno/Release
echo "DONE copy GenomicsDB executables to ${host_dest_dir} .... container name is ${cntr_name}"