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
#
# The code utilized some features of RedHat pre-built container
FROM centos:7  

LABEL vendor="Intel Corporation" name="VCF Combiner" description="Combine VCF files into a single VCF file"

RUN yum install -y --setopt=tsflags=nodocs epel-release &&  \
    yum repolist && \
    yum install -y python34.x86_64 && \
    yum install -y python34-setuptools && \
    easy_install-3.4 pip && \
    yum install -y --setopt=tsflags=nodocs libcsv mpich openssl zlib unzip.x86_64 && \  
    yum clean all && \
    pip3 install PyVCF

ADD ./usr /usr 
ADD ./etc /etc
ADD ./root /root
COPY tmp/vcf2tiledb /usr/bin/

WORKDIR /tmp

ENV BASH_ENV=/etc/profile.d/cont-env.sh HOME=/home/default PATH=$PATH:/usr/lib64/mpich/bin 

RUN  groupadd -r default -f -g 5658 && \
     useradd -u 5658 -g default -o -c "Default User" default -s /sbin/nologin

#USER default

ENTRYPOINT ["/usr/bin/container-entrypoint"]

CMD ["container-usage"]
