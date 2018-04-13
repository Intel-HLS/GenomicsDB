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
#
# The code utilized some features of RedHat pre-built container
[[ $# -eq 1 ]] && cache_opt='--no-cache' || cache_opt=''

PRODUCT_LABEL="name=\"$PRODUCT_NAME\" version=\"$PRODUCT_VERSION\" description=\"$PRODUCT_DESCRIPTION\""
which $TARGET_EXEC >/dev/null 2>&1 || echo "Could not find vcf2tiledb", exit 1
echo "$TARGET_EXEC version is $(${TARGET_EXEC} --version)"
full_version="$(${TARGET_EXEC} --version)-$PRODUCT_VERSION"
echo "START building $PRODUCT_NAME version $full_version docker..."

DOCKER_PNAME=`python -c "x='${PRODUCT_NAME}'.lower().split();print('%s_%s' % (x[0], x[1]))"`
pushd $DOCKER_ROOT >/dev/null 2>&1
[[ -d $DOCKER_PNAME ]] && rm -r $DOCKER_PNAME
mkdir -p $DOCKER_PNAME/query_interval_samples
cp -r $DOCKER_GRAND_ROOT/geno_docker_template/* $DOCKER_PNAME/
cp $DOCKER_GRAND_ROOT/known_vid_mapper_files/* $DOCKER_PNAME/usr/share/cont-intel/
# copy test_resource
mkdir -p $DOCKER_PNAME/usr/share/cont-intel/test_res
find $(dirname $MY_TEST_RESOURCE) -maxdepth 1 -type f -exec cp -t $DOCKER_PNAME/usr/share/cont-intel/test_res {} +
cp $MY_TEST_RESOURCE/* $DOCKER_PNAME/usr/share/cont-intel/test_res

[[ -f Dockerfile ]] && cp Dockerfile $DOCKER_PNAME/ || cp $DOCKER_GRAND_ROOT/Dockerfile $DOCKER_PNAME/
sed -i "s/<label>/$PRODUCT_LABEL/" $DOCKER_PNAME/Dockerfile

cp *.py $DOCKER_PNAME/usr/bin/
if [ -f "$DOCKER_PNAME/usr/bin/$DOCKER_PNAME.py" ]
then
    pushd "$DOCKER_PNAME/usr/bin"
    ln -s "$DOCKER_PNAME.py" $DOCKER_PNAME
    popd
fi
chmod +x -R $DOCKER_PNAME/usr
mkdir -p $DOCKER_PNAME/tmp/bin
cp $(which $TARGET_EXEC) $DOCKER_PNAME/tmp/bin/
genolib=$(ls -t -1 $(dirname ${DOCKER_ROOT})/PyGenomicsDBPerf-*-py3-none-any.whl | head -n 1)
cp $genolib $DOCKER_PNAME/tmp/
sed -i "s|<genolib>|$(basename ${genolib})|g" $DOCKER_PNAME/Dockerfile
[[ -f help.txt ]] && cp help.txt $DOCKER_PNAME/usr/share/cont-docs/

cd $DOCKER_PNAME
docker build -t $DOCKER_PNAME:$full_version '--no-cache' --build-arg http_proxy=$http_proxy --build-arg https_proxy=$http_proxy --build-arg HTTP_PROXY=$http_proxy --build-arg HTTPS_PROXY=$http_proxy .
[[ $? -eq 0 ]] && echo "INFO: created docker image $DOCKER_PNAME:$full_version" || echo "INFO: failed to create docker image $DOCKER_PNAME:$full_version"
popd >/dev/null 2>&1
