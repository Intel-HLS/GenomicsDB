#! /bin/bash

# -----------------------------------------------------------------------------

# The MIT License (MIT)
# Copyright (c) 2016 Intel Corporation
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

# -----------------------------------------------------------------------------

create_load_tiledb_cfg() {
    outfile="${INITDIR}/load_to_tiledb.cfg"
    rm ${outfile}
    touch ${outfile}
    if [ 0 ]
    then
        echo "[loader]"
        echo "num_processes = 1"
        echo "executable = ${GDBDIR}/bin/vcf2tiledb"
        echo "tile_loader_json = ${INITDIR}/tile_loader.json"
    fi >> ${outfile}
}

# -----------------------------------------------------------------------------

create_work_space() {
    rm -rf "${INITDIR}"
    cp -r ${OLD_INITDIR} ${INITDIR}
    create_tiledb_workspace "${INITDIR}/workspace"
    mkdir -p "${INITDIR}/outputs"
    (
        cd ${INITDIR}/inputs/vcfs
        tabix -f t0_1_2_combined.vcf.gz
    )
}

# -----------------------------------------------------------------------------

create_vcf_import_cfg() {
    outfile="${INITDIR}/vcf_import.config"
    rm ${outfile}
    touch ${outfile}
    if [ 0 ]
    then
        echo "{"
        echo "\"workspace\": \"${INITDIR}/workspace\","
        echo "\"vcf_type\": \"composite\","
        echo "\"array\": \"sql_mapper_test\","
        echo "\"assembly\": \"hg19\","
        echo "\"dburi\": \"postgresql+psycopg2://postgres:postgres@localhost:5432/gendb\""
        echo "}"
    fi >> ${outfile}
}

# -----------------------------------------------------------------------------

vcf_only_metadata() {
    echo "python ${SCRIPT} -c ${INITDIR}/vcf_import.config -d ${INITDIR}/outputs -i ${INITDIR}/inputs/vcfs/t0_1_2_combined.vcf.gz -a ${INITDIR}/inputs/callsets/t0_1_2.json"
    echo -n "Press any key to continue: "
    read prompt
    python ${SCRIPT} -c ${INITDIR}/vcf_import.config -d ${INITDIR}/outputs -i ${INITDIR}/inputs/vcfs/t0_1_2_combined.vcf.gz -a ${INITDIR}/inputs/callsets/t0_1_2.json
}

# -----------------------------------------------------------------------------

#vcf_metadata_and_tile() {
    #echo "python ${SCRIPT} -c ${INITDIR}/vcf_import.config -d ${INITDIR}/outputs -i ${INITDIR}/inputs/vcfs/t0_1_2_combined.vcf.gz -a ${INITDIR}/inputs/callsets/t0_1_2.json -l ${INITDIR}/load_to_tile.cfg"
    #echo -n "Press any key to continue: "
    #read prompt
    #python ${SCRIPT} -c ${INITDIR}/vcf_import.config -d ${INITDIR}/outputs -i ${INITDIR}/inputs/vcfs/t0_1_2_combined.vcf.gz -a ${INITDIR}/inputs/callsets/t0_1_2.json -l ${INITDIR}/load_to_tile.cfg
#}

# -----------------------------------------------------------------------------

vcf_only_tile() {
    echo "python ${INITDIR}/run.py ${GDBDIR} ${GDBDIR}"
    #echo "COMMENTED OUT THE ABOVE COMMAND DUE TO VcfAdapterException - VERIFY"
    python ${INITDIR}/run.py ${GDBDIR} ${GDBDIR}
}

# -----------------------------------------------------------------------------

init_venv() {
    echo "Ensure all prerequiste packages are installed."
    echo "For e.g., libpq-dev python-dev postgresql-server-dev etc"
    echo -n "Are all the prerequisites installed (y/n)?: "

    read answer
    if [ "${answer}" != "y" ]
    then
        echo "Please install all prerequisites. Exiting"
        echo "--------------------------------------"
        exit
    fi

    echo "--------------------------------------"
    cd ${HOME}
    virtualenv venv
    source venv/bin/activate

    export PATH="${GDBDIR}/bin:${PATH}"
    cd ${APIDIR}
    pip install -r requirements.txt
    python setup.py develop
    export PYTHONPATH="${APIDIR}:${APIDIR}/web"

    echo "--------------------------------------"
    echo "Using postgresql DB user postgres. Enter password when requested"
    createdb -U postgres gendb
    createlang -U postgres plpgsql gendb
    cd metadb

    if [ ! -f "alembic.ini" ]
    then
        echo "Make sure alembic.ini is copied to `pwd` and edited appropriately"
        exit
    fi
}

# -----------------------------------------------------------------------------

if [ -z ${GDBDIR}  -o  -z ${APIDIR} ]
then
    echo "GDBDIR and APIDIR must be defined"
    echo "Usage: <load_genomics_metadata.sh GenomicsDB GenomicsSampleAPIs>"
    exit
fi

if [ -z ${http_proxy}  -o  -z ${https_proxy} ]
then
    echo "http_proxy and https_proxy must be defined"
    exit
fi

OLD_INITDIR="${GDBDIR}/tests"
INITDIR="/tmp/sql_mapper"
SCRIPT="${APIDIR}/utils/vcf2tile.py"

echo "--------------------------------------"
echo "RSM: SHELL - GDBDIR - <${GDBDIR}>"
echo "RSM: SHELL - APIDIR - <${APIDIR}>"
echo "RSM: SHELL - INITDIR - <${INITDIR}>"
echo "RSM: SHELL - SCRIPT - <${SCRIPT}>"
echo "--------------------------------------"

init_venv
echo -n "Is Schema Initialized (y/n)?: "
read answer
if [ "${answer}" == "n" ]
then
    # load schema into DB
    #alembic init alembic
    #alembic revision -m "create field_set table"
    alembic upgrade head

    error=$?
    if [ $error -ne 0 ]
    then
        echo "ERROR: <${error}>"
        exit
    fi

    create_work_space
    create_load_tiledb_cfg
    create_vcf_import_cfg
fi

vcf_only_metadata
#vcf_metadata_and_tile
#vcf_only_tile

echo "--------------------------------------"
# -----------------------------------------------------------------------------

