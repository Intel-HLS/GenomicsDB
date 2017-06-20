### GenomicsDB Building Docker

The centos-build is the source code for building a GenomicsDB building environment Docker image. The details about GenomicsDB building environment can be found at [GenomicsDB Wiki]( https://github.com/Intel-HLS/GenomicsDB/wiki/Compiling-GenomicsDB#building.)

To create a GenomicsDB building environment Docker image, run the following command:

<code>docker build -t name_of_the_builder_image --no-cache .</code>

Once the image is created. You can build the latest GenomicsDB executables by running:

<code>docker run -it -v /path/2/output/:/output/ [-e  GDB_BRANCH=branch_name_default_is_master] name_of_the_builder_image</code>
