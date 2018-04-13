#! /usr/bin/python3

"""
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
"""

import os
import os.path

MPIRUN_PATH = "/usr/lib64/mpich/bin/mpirun"
DEFAULT_DB_NAME = "genomicsdb_loader.db"
DEFAULT_WORKSPACE = os.path.expanduser("~")
NO_VERSION = 'not available'

L_MARKER = "LLL"
Q_MARKER = "QQQ"

LOGGER_ROOT = os.path.join(os.environ['HOME'], 'GenomicsDBPerfTest', 'run_logs', 't_run_logs', 'logfiles')

GDB_COMMANDS = ['vcf2tiledb', 'gt_mpi_gather', 'java']
LOADER_EXEC = GDB_COMMANDS[0]
CHECK_VERSION_LIST = GDB_COMMANDS[:2]

DENSITY_NAMES = ['dense', 'meanmid', 'sparse']
DevNull = open(os.devnull, 'wb', 0)

WAIT_TIME_FOR_LAUNCHING = 0.5           # 0.5 sec
DEFAULT_SAMPLING_INTERVAL = 1                    # in sec

CONFIG_DIRNAME = 'config_files'
RUN_SCRIPT = "main_remote.py"

ENV_TILEDB_ROOT = 'TILEDB_ROOT'  #?
DEFAULT_TDB_PREFIX = "/mnt/app_hdd1/scratch/mingperf/tiledb-ws"
GT_MPI_GATHER_OPTIONS = ['--print-calls', '--produce-Broad-GVCF', '--print-AC', '--produce-interesting-positions']

GENOMICSDB_TIMER = 'GENOMICSDB_TIMER'

ENV_LIB_PATH = 'LD_LIBRARY_PATH'
DEFAULT_LIB_PATHS = [
    os.path.join(os.environ['HOME'], 'opt', 'protobuf', 'lib'),
    '/usr/lib64/mpich/lib/',
    '/opt/rh/devtoolset-4/root/usr/lib64',
    '/opt/rh/devtoolset-4/root/usr/lib']

CLEAN_CACHE_SCRIPT_PATH = "/tools/utilities/bin/clear_caches.sh"

TIME_FORMATTER = "0~%C,1~%e,2~%P,3~%F,4~%R,5~%w,6~%I,7~%O,8~%c,9~%x,10~%M,11~%t,12~%K,13~%U,14~%S,15~%D,16~%Z"

TIME_ROW_LABELS = [
    'Command',
    'Elapsed (wall clock) time (sec)',
    'Percent of CPU this job got %',
    'Major (requiring I/O) page faults',
    'Minor (reclaiming a frame) page faults',
    'Involunteer context switches',
    'File system input',
    'File system output',
    'Volunteer context switches',
    'Exit status',
    'Max resident set size (kbytes)',
    'Avg resident set size (kbytes)',
    'Average total size (kbytes)',
    'User CPU-sec time',
    'System CPU-sec time',
    'Average unshared data size (kbytes)',
    'Page size (bytes)']

# string spilled by vcf2tile and gt_mpi_gather
genome_profile_tags = {'fetch from vcf': 'fv', 'combining cells':'cc', 'flush output': 'fo',\
    'sections time': 'st', 'time in single thread phase()': 'ts', 'time in read_all()': 'tr'}
genome_queryprof_tags = {'genomicsdb cell fill timer': 'cf', 'bcf_t serialization': 'bs', \
   'operator time': 'ot', 'sweep at query begin position': 'sq', 'tiledb iterator': 'ti', \
   'total scan_and_produce_broad_gvcf time for rank 0': 'tt', 'tiledb to buffer cell': 'tb', \
   'bcf_t creation time': 'bc'}
