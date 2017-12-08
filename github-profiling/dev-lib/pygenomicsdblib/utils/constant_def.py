#! /usr/bin/python3
import os
import os.path
import platform
from collections import namedtuple

MPIRUN_PATH = "/usr/lib64/mpich/bin/mpirun"
DEFAULT_DB_NAME = "genomicsdb_loader.db"
DEFAULT_WORKSPACE = os.path.expanduser("~")
NO_VERSION = 'not available' 

L_MARKER = "LLL" 
Q_MARKER = "QQQ"

LOGGER_ROOT = os.path.join(os.environ['HOME'], 'GenomicsDBPerfTest', 'run_logs', 't_run_logs', 'logfiles') if platform.system() == 'Linux' else '/path/to/run_logs'

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
DEFAULT_TDB_PREFIX = "/path/to/tiledb-ws"
GT_MPI_GATHER_OPTIONS = ['--print-calls', '--produce-Broad-GVCF', '--print-AC', '--produce-interesting-positions']

GENOMICSDB_TIMER = 'GENOMICSDB_TIMER'

ENV_LIB_PATH = 'LD_LIBRARY_PATH'
DEFAULT_LIB_PATHS = [
    os.path.join(os.environ['HOME'], 'opt', 'protobuf', 'lib',
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

#### config json templates
# 'to_tiledb' is based on 1000g
# "size_per_column_partition": num_vcfs * COL_PARTITION_SIZE_UNIT, for 1000g is 1000
__COL_PARTITION_SIZE_UNIT = 16384                     # used in load_cfg:size_per_column_partition
__G1000_COL_PARTITION_SIZE_UNIT = 1000 * __COL_PARTITION_SIZE_UNIT  # genomics 1000
__GENODB_ATTRIBUTES = ["REF", "ALT", "BaseQRankSum", "MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", \
"GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP"]
__ILLMN_GENODB_ATTRIBUTES = ["REF", "ALT", "QUAL", "GT"]
__cfg_json_t = { 
    "combine_vcf" : {
        "produce_combined_vcf": True,
        "num_parallel_vcf_files": 1,
        "callset_mapping_file" : "",
        "vcf_output_format" : "z",
        "vcf_output_filename" : "",
        "reference_genome" : "", 
        "column_partitions" : [],
        "do_ping_pong_buffering" : True,
        "offload_vcf_output_processing" : True,
        "produce_GT_field" : True,
        "size_per_column_partition" : __COL_PARTITION_SIZE_UNIT,
        "vid_mapping_file": ""},
    "to_tiledb" :  {
        "callset_mapping_file": "",
        "column_partitions": [],
        "compress_tiledb_array": True,
        "delete_and_create_tiledb_array": True,
        "disable_synced_writes": True,
        "discard_vcf_index": True,
        "do_ping_pong_buffering": False,
        "ignore_cells_not_in_partition": False,
        "num_cells_per_tile": 1000,
        "num_parallel_vcf_files": 1,
        "offload_vcf_output_processing": True,
        "produce_combined_vcf": False,
        "produce_tiledb_array": True,
        "reference_genome": "/path/to/Homo_sapiens_assembly19.fasta",
        "row_based_partitioning": False,
        "segment_size": 10485760,
        "size_per_column_partition": __G1000_COL_PARTITION_SIZE_UNIT,
        "tiledb_compression_level": 6,
        "treat_deletions_as_intervals": False,
        "lb_callset_row_idx": 0,
        "ub_callset_row_idx": 999,
        "vcf_header_filename": "",
        "vcf_output_format": "z",
        "vid_mapping_file": "/path/to/vid.json"},
    "query_full_scan" : {
        "array": ["TEST0"],
        "callset_mapping_file": "temp",
        "query_attributes": __GENODB_ATTRIBUTES,
        "reference_genome": "/path/to/Homo_sapiens_assembly19.fasta",
        "scan_full": True,
        "segment_size": 10000,
        "vid_mapping_file": "/path/to/vid.json",
        "workspace": ""},
    "query_density" : {
        "array": ["TEST0"],
        "callset_mapping_file": "temp",
        "query_attributes": __GENODB_ATTRIBUTES,
        "reference_genome": "/path/to/broad_reference/Homo_sapiens_assembly19.fasta",
        "query_column_ranges" : [[3095999, 9287998, 15479997]],
        "segment_size": 10000,
        "vid_mapping_file": "/path/to/vid.json",
        "workspace": ""},
    "histogram" : {
        "column_partitions" : [{"begin": 0}],
        "callset_mapping_file" : "",
        "vid_mapping_file" : "",
        "num_converter_processes" : 1,
        "max_histogram_range" : 3999999999,
        "num_bins" : 100000,
        "size_per_column_partition": 16384,
        "num_parallel_vcf_files" : 1}
}
def __get_nt(name, cfg):
    nt_tags = namedtuple(name, ','.join(cfg.keys()))
    return nt_tags(**cfg)

def __query_nt(name, qcfg, num_parallel):
    qcfg["array"] = ["TEST%d" % i for i in range(num_parallel)]
    return __get_nt(name, qcfg)
    
LOADER_CFG_T = lambda: __get_nt('nt_tdb_loader', __cfg_json_t['to_tiledb']) 
# partition is empty
LOADER_COMBINER_CFG_T = lambda: __get_nt('nt_combine_vcf', __cfg_json_t['combine_vcf'])
# query_column_ranges is empty
QUERY_CFG_T = lambda num_par: __query_nt("nt_query", __cfg_json_t["query_density"], num_par)
QUERY_FULL_SCAN_CFG_T = lambda num_par: __query_nt("nt_query_full", __cfg_json_t["query_full_scan"], num_par)

CREATE_COL_PARTITION_F = lambda array_start_pos, tdb_ws: [{"array": "TEST%d" % (i), "begin": begin, "workspace": tdb_ws} \
            for i, begin in enumerate(array_start_pos)]

HISTOGRAM_CFG_T = lambda: __get_nt('nt_histogram', __cfg_json_t['histogram'])
