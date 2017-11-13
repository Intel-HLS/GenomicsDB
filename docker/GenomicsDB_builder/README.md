### GenomicsDB Building Docker

The MIT License (MIT)
Copyright (c) 2016-2017 Intel Corporation

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

--------------------

Building GenomicsDB could be time-consuming. This project create a docker image that 
builds environment for building Intel GenomicsDB. Running the docker image builds GenomicsDB 
utilities. The details about building GenomicsDB can be found at [GenomicsDB Wiki]( https://github.com/Intel-HLS/GenomicsDB/wiki/Compiling-GenomicsDB#building.)

The centos:7 is the base image of GenomicsDB building environment Docker image.

##### Files:
* dp_builder.bash - shell script for building the docker image
* dp_run_builder - shell script for running the built docker image
* cont-gdb_builder.bash - shell script for building GenomicsDB
* enabledevtoolset-4.sh - shell script for enable centos dev tools.
* GDB-builder_centos/Dockerfile - docker script
* GDB-builder_centos/other files are part of RedHat template.

##### Usage:

1. set an environment variable MY_DOCKER_ROOT that point to your_GenomicsDB_root/docker/GenomicsDB_builder
2. run <b>'dp_builder.bash'</b> this script creates a docker image name 'genomicsdb_builder'. The docker image sets up environment for building GenomicsDB.
3. set an environment variable GENOME_VOLUME=genome-shared. Docker 'genomicsdb_builder' uses docker volume $GENOME_VOLUME, where it stores the outputs.   
4. run <b>'dp_run_builder.bash [branch name]'</b>, 
    where the 'branch name' is a GenommicsDB github branch. The default branch name is 'master'
    this script pulls the latest version of the selected branch from github, then builds the images.
5. Once the image is built, just run "dp_run_builder.bash [branch name]" to build GenomicsDB. The instances are named for easy backtracking.
