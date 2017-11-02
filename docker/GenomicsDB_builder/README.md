### GenomicsDB Building Docker
<code>
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
</code>
Building GenomicsDB could be time-consuming. This project create a docker image that 
builds environment for building Intel GenomicsDB. Running the docker image builds GenomicsDB 
utilities. The details about building GenomicsDB can be found at [GenomicsDB Wiki]( https://github.com/Intel-HLS/GenomicsDB/wiki/Compiling-GenomicsDB#building.)

The centos:7 is the base image of GenomicsDB building environment Docker image.

The 2 bash scripts are used for creating docker image and running instance respectively. The scripts are suitable for our environment. You may need adapt the scripts to your environment, especially the network part.

Usage:

1. <code>run 'dp_builder.bash'</code>   
    this script creates a docker image name 'genomicsdb_builder'
2. <code>run 'dp_run_builder.bash [branch name]', 
    where the 'branch name' is a GenommicsDB github branch. The default branch name is 'master'</code>
     this script runs the docker image 'genomicsdb_builder' that builds GenomicdDB software. The docker instance is name as "GDB-{branch_name}_{timestamp}", example "GDB-branch1_171030_1230"
3. the rest code are part of docker building.

The docker instance normally are named with generated UUID. User named instances make identify the instances easy, which is useful for backtracking.
