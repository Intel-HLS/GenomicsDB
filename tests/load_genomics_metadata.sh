#! /bin/bash
# Ramesh Mantri
# Thu Apr 13 15:33:02 PDT 2017

# -----------------------------------------------------------------------------

create_load_tile_cfg() {
    outfile="${INITDIR}/load_to_tile.cfg"
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
    rm -rf "${INITDIR}/workspace"
    #mkdir -p "${INITDIR}/workspace"
    create_tiledb_workspace "${INITDIR}/workspace"
    mkdir -p "${INITDIR}/outputs"
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
        echo "\"array\": \"rsm_test\","
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

vcf_metadata_and_tile() {
    echo "python ${SCRIPT} -c ${INITDIR}/vcf_import.config -d ${INITDIR}/outputs -i ${INITDIR}/inputs/vcfs/t0_1_2_combined.vcf.gz -a ${INITDIR}/inputs/callsets/t0_1_2.json -l ${INITDIR}/load_to_tile.cfg"
    echo -n "Press any key to continue: "
    read prompt
    python ${SCRIPT} -c ${INITDIR}/vcf_import.config -d ${INITDIR}/outputs -i ${INITDIR}/inputs/vcfs/t0_1_2_combined.vcf.gz -a ${INITDIR}/inputs/callsets/t0_1_2.json -l ${INITDIR}/load_to_tile.cfg
}

# -----------------------------------------------------------------------------

vcf_only_tile() {
    echo "python ${INITDIR}/run.py ${GDBDIR} ${GDBDIR}"
    #echo "COMMENTED OUT THE ABOVE COMMAND DUE TO VcfAdapterException - VERIFY"
    python ${INITDIR}/run.py ${GDBDIR} ${GDBDIR}
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

INITDIR="${GDBDIR}/tests"
SCRIPT="${APIDIR}/utils/vcf2tile.py"

echo "--------------------------------------"
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

echo "RSM: SHELL - GDBDIR - <${GDBDIR}>"
echo "RSM: SHELL - APIDIR - <${APIDIR}>"
echo "RSM: SHELL - INITDIR - <${INITDIR}>"
echo "RSM: SHELL - SCRIPT - <${SCRIPT}>"

echo "--------------------------------------"
echo "Using postgresql DB user postgres. Enter password when requested"
createdb -U postgres gendb
createlang -U postgres plpgsql gendb
cd metadb

if [ ! -f "alembic.ini" ]
then
    echo "Make sure alembic.ini is copied to `pwd` and edited appropriately:"
    exit
fi

alembic upgrade head
create_work_space
create_load_tile_cfg
create_vcf_import_cfg

echo "--------------------------------------"
vcf_only_metadata
#vcf_metadata_and_tile
vcf_only_tile
echo "--------------------------------------"
# -----------------------------------------------------------------------------

