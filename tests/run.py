#!/usr/bin/env python

import json
import tempfile
import subprocess
import hashlib
import os
import sys
import shutil

query_json_template_string="""
{   
        "workspace" : "",
        "array" : "",
        "vcf_header_filename" : ["inputs/template_vcf_header.vcf"],
        "query_column_ranges" : [ [ [0, 10000000000 ] ] ],
        "query_row_ranges" : [ [ [0, 2 ] ] ],
        "reference_genome" : "inputs/chr1_10MB.fasta.gz",
        "query_attributes" : [ "REF", "ALT", "BaseQRankSum", "MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP" ]
}"""

def create_query_json(ws_dir, test_name, query_column_range):
    test_dict=json.loads(query_json_template_string);
    test_dict["workspace"] = ws_dir
    test_dict["array"] = test_name
    test_dict["query_column_ranges"] = [ [ query_column_range ] ]
    return test_dict;


loader_json_template_string="""
{
    "row_based_partitioning" : false,
    "column_partitions" : [
        {"begin": 0, "workspace":"", "array": "" }
    ],
    "callset_mapping_file" : "",
    "vid_mapping_file" : "inputs/vid.json",
    "size_per_column_partition": 3000,
    "treat_deletions_as_intervals" : true,
    "vcf_header_filename": "inputs/template_vcf_header.vcf",
    "reference_genome" : "inputs/chr1_10MB.fasta.gz",
    "num_parallel_vcf_files" : 1,
    "do_ping_pong_buffering" : false,
    "offload_vcf_output_processing" : false,
    "discard_vcf_index": true,
    "produce_combined_vcf": true,
    "produce_tiledb_array" : true,
    "delete_and_create_tiledb_array" : true,
    "compress_tiledb_array" : false,
    "segment_size" : 1048576,
    "num_cells_per_tile" : 10
}""";

def create_loader_json(ws_dir, test_name):
    test_dict=json.loads(loader_json_template_string);
    test_dict["column_partitions"][0]["workspace"] = ws_dir;
    test_dict["column_partitions"][0]["array"] = test_name;
    test_dict["callset_mapping_file"] = 'inputs/callsets/'+test_name+'.json';
    return test_dict;

def cleanup_and_exit(tmpdir, exit_code):
    shutil.rmtree(tmpdir, ignore_errors=True)
    sys.exit(exit_code);

def main():
    #Switch to tests directory
    parent_dir=os.path.dirname(os.path.realpath(__file__))
    os.chdir(parent_dir)
    #Zero line coverage
    subprocess.call('lcov --directory ../ --zerocounters', shell=True);
    exe_path = '../bin/'
    tmpdir = tempfile.mkdtemp()
    ws_dir=tmpdir+os.path.sep+'ws';
    loader_tests = [
            { "name" : "t0_1_2", 'stdout_md5sum_hash' : '6063db9832fcdeb25d24dcc0630ec499',
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000], "stdout_md5sum_hash": [
                        "0e13f0030b56e906eda9124e3c922949",
                        "976e1b4cb258bf791206e60e890c0a48",
                        "3b9b56b9a9a64199fe20d919703aca98"
                    ] },
                    { "query_column_ranges" : [12150, 1000000000], "stdout_md5sum_hash": [
                        "fc7dee1c356ed054b8858bbe7b8da83e",
                        "976e1b4cb258bf791206e60e890c0a48",
                        "76d636ec426fc7cf3b53508554a95224"
                        ] }
                    ]
            },
            { "name" : "t0_overlapping" },
            { "name" : "t6_7_8", 'stdout_md5sum_hash' : '6ddaed1219aae6ac2f4ea5da3f312178',
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000], "stdout_md5sum_hash": [
                        "a8e66cc7df0001b64da650d3a901d6ef",
                        "b677b16582b477cef2fb3b8e71417f5b",
                        "6ddaed1219aae6ac2f4ea5da3f312178"
                    ] },
                    { "query_column_ranges" : [8029500, 1000000000], "stdout_md5sum_hash": [
                        "81427830d0d43701b49b6c22bd9e7369",
                        "b677b16582b477cef2fb3b8e71417f5b",
                        "19d86392bf34a3f737bb47a99cb1b8d5"
                    ] }
                    ]
            },
    ];
    for test_params_dict in loader_tests:
        test_name = test_params_dict['name']
        test_loader_dict = create_loader_json(ws_dir, test_name);
        if(test_name == "t0_overlapping"):
            test_loader_dict["produce_combined_vcf"] = False;
        loader_json_filename = tmpdir+os.path.sep+test_name+'.json'
        with open(loader_json_filename, 'wb') as fptr:
            json.dump(test_loader_dict, fptr);
            fptr.close();
        pid = subprocess.Popen(exe_path+os.path.sep+'vcf2tiledb '+loader_json_filename, shell=True,
                stdout=subprocess.PIPE);
        stdout_string = pid.communicate()[0]
        if(pid.returncode != 0):
            sys.stderr.write('Loader test: '+test_name+'failed\n');
            cleanup_and_exit(tmpdir, -1);
        md5sum_hash_str = str(hashlib.md5(stdout_string).hexdigest())
        if('stdout_md5sum_hash' in test_params_dict and md5sum_hash_str != test_params_dict['stdout_md5sum_hash']):
            sys.stderr.write('Loader stdout mismatch for test: '+test_name+'\n');
            cleanup_and_exit(tmpdir, -1);
        if('query_params' in test_params_dict):
            for query_param_dict in test_params_dict['query_params']:
                test_query_dict = create_query_json(ws_dir, test_name, query_param_dict["query_column_ranges"])
                query_types_list = [ ('calls','--print-calls'), ('variants',''), ('vcf','--produce-Broad-GVCF') ]
                idx = 0;
                for query_type,cmd_line_param in query_types_list:
                    query_json_filename = tmpdir+os.path.sep+test_name+'_'+query_type+'.json'
                    with open(query_json_filename, 'wb') as fptr:
                        json.dump(test_query_dict, fptr);
                        fptr.close();
                    pid = subprocess.Popen(exe_path+os.path.sep+'gt_mpi_gather -l '+loader_json_filename+' -j '
                            +query_json_filename+' '+cmd_line_param, shell=True,
                            stdout=subprocess.PIPE);
                    stdout_string = pid.communicate()[0]
                    if(pid.returncode != 0):
                        sys.stderr.write('Query test: '+test_name+'-'+query_type+' failed\n');
                        cleanup_and_exit(tmpdir, -1);
                    md5sum_hash_str = str(hashlib.md5(stdout_string).hexdigest())
                    if('stdout_md5sum_hash' in query_param_dict and md5sum_hash_str !=
                            query_param_dict['stdout_md5sum_hash'][idx]):
                        sys.stderr.write('Mismatch in query test: '+test_name+'-'+query_type+'\n');
                        cleanup_and_exit(tmpdir, -1);
                    idx += 1
    coverage_file='coverage.info'
    subprocess.call('lcov --directory ../ --capture --output-file '+coverage_file, shell=True);
    subprocess.call("lcov --remove "+coverage_file+" '/opt*' '/usr*' 'dependencies*' -o "+coverage_file, shell=True);
    cleanup_and_exit(tmpdir, 0); 

if __name__ == '__main__':
    main()
