### GenomicsDB Building Docker

Building GenomicsDB could be time-consuming. The project uses docker container for building GenomicsDB. We share the project as an example of building Intel [GenomicsDB]( https://github.com/Intel-HLS/GenomicsDB/wiki/Compiling-GenomicsDB#building.) project. The project was tested for building GenomicsDB 'master' branch on Nov 16, 2017.

This project provides 2 bash files and associated docker files. 'dp_builder.bash' creates a docker image that builds environment for building Intel GenomicsDB. 'dp_run_builder' builds GenomicsDB utilities with the current version of the branch (the default is 'master') pulled from [GenomicsDB github](https://github.com/Intel-HLS/GenomicsDB).

<b>Adapt the scripts accordingly</b> to your needs and GenomicsDB changes. The details about building GenomicsDB can be found at [GenomicsDB Wiki]( https://github.com/Intel-HLS/GenomicsDB/wiki/Compiling-GenomicsDB#building.)

The centos:7 is the base image of GenomicsDB building environment Docker image.

##### Usage:
1. set an environment variable GENOME_VOLUME, such as 'GENOME_VOLUME=genome-shared'.
2. optional, set an environment variable MY_DOCKER_ROOT that point to your_GenomicsDB_root/docker/GenomicsDB_builder folder. If MY_DOCKER_ROOT is not set, it assumes that you are in your_GenomicsDB_root/docker/GenomicsDB_builder folder. 
3. run <b>'dp_builder.bash'</b>:
    * the script creates a docker image name 'genomicsdb_builder'. 
    * once the image is built, just run "dp_run_builder.bash [branch name]" to build GenomicsDB.
    * need to re-build the image if the Dockerfile is changed.
4. run <b>'dp_run_builder.bash [branch_name]'</b>:
    * where the 'branch name' is a name of a GenommicsDB github branch. The default branch name is 'master'.
    * the script launches a named docker instance of docker image 'genomicsdb_builder'. The name of the docker instance contains branch name and timestamp, such as 'GDB-master_171116_1830'. 
    * the docker instance pulls the latest version of the selected branch from github, then builds the GenommicsDB libraries and utilities.
    * then, the script copies the build to ${HOME}/docker-geno_[branch_name]_[builddate] (for example: ${HOME}/docker-geno_columnar_field_171117)

##### Files you may want or need modify:
* dp_builder.bash - shell script for building the docker image.
* dp_run_builder - shell script for running the built docker image.
* build_gdb.bash - the script that builds GenomicsDB. <b>Modify line with cmake</b> to your needs according to instruction described at [GenomicsDB Wiki]( https://github.com/Intel-HLS/GenomicsDB/wiki/Compiling-GenomicsDB#building.)
* GDB-builder_centos/Dockerfile - docker script. <b>Add/remove yum packages</b> in this file when needed.

##### Files are more static:
* enabledevtoolset-4.sh - shell script for enable centos dev tools.
* GDB-builder_centos/other files are part of RedHat template.

---------------------------

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
