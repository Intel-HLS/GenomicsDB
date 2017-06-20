### GenomicsDB Building Docker

The centos-build is the source code for building a GenomicsDB building environment Docker image. The details about GenomicsDB building environment can be found at [GenomicsDB Wiki]( https://github.com/Intel-HLS/GenomicsDB/wiki/Compiling-GenomicsDB#building.)

To create a GenomicsDB building environment Docker image, run the following command:

<code>docker build -t name_of_the_builder_image --no-cache .</code>

Once the image is created. You can build the latest GenomicsDB executables by running:

<code>docker run -it -v /path/2/output/:/output/ [-e  GDB_BRANCH=branch_name_default_is_master] name_of_the_builder_image</code>

In this case, the docker container builds GenomicsDB with the following options:

* CMAKE_BUILD_TYPE=Release
* DO_PROFILING=False
* DISABLE_OPENMP=True
* BUILD_JAVA=False
* DO_PROFILING=False  
* CMAKE_INSTALL_PREFIX=/output/
* PROTOBUF_LIBRARY=dont_worry_the_docker_manage_it

You can pass GenomicsDB building options via docker run command arguments. The docker will pass your arguments to cmake. Since the docker builds protobuf library internally, you just ignore the PROTOBUF_LIBRARY option.
