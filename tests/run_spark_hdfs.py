#!/usr/bin/env python

#The MIT License (MIT)
#Copyright (c) 2018 University of California, Los Angeles and Intel Corporation

#Permission is hereby granted, free of charge, to any person obtaining a copy of 
#this software and associated documentation files (the "Software"), to deal in 
#the Software without restriction, including without limitation the rights to 
#use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
#the Software, and to permit persons to whom the Software is furnished to do so, 
#subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all 
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
#FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
#COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
#IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
#CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import json
import tempfile
import subprocess
import hashlib
import os
import sys
import shutil
import difflib
from collections import OrderedDict

query_json_template_string="""
{   
        "workspace" : "",
        "array" : "",
        "vcf_header_filename" : ["inputs/template_vcf_header.vcf"],
        "query_column_ranges" : [ [ [0, 10000000000 ] ] ],
        "query_row_ranges" : [ [ [0, 3 ] ] ],
        "query_block_size" : 10000,
        "query_block_size_margin" : 500,
        "reference_genome" : "inputs/chr1_10MB.fasta.gz",
        "query_attributes" : [ "REF", "ALT", "BaseQRankSum", "MQ", "RAW_MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP", "PID", "PGT" ]
}"""

vcf_query_attributes_order = [ "END", "REF", "ALT", "BaseQRankSum", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "MQ", "RAW_MQ", "MQ0", "DP", "GT", "GQ", "SB", "AD", "PL", "PGT", "PID", "MIN_DP", "DP_FORMAT", "FILTER" ];
query_attributes_with_DS_ID = [ "REF", "ALT", "BaseQRankSum", "MQ", "RAW_MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP", "PID", "PGT", "DS", "ID" ];
query_attributes_with_PL_only = [ "PL" ]
query_attributes_with_MLEAC_only = [ "MLEAC" ]
default_segment_size = 40

def create_query_json(ws_dir, test_name, query_param_dict, test_dir):
    test_dict=json.loads(query_json_template_string);
    test_dict["workspace"] = ws_dir
    test_dict["array"] = test_name
    test_dict["query_column_ranges"] = [ [ query_param_dict["query_column_ranges"] ] ]
    if("vid_mapping_file" in query_param_dict):
        test_dict["vid_mapping_file"] = query_param_dict["vid_mapping_file"];
    if("callset_mapping_file" in query_param_dict):
        test_dict["callset_mapping_file"] = query_param_dict["callset_mapping_file"];
    if("query_attributes" in query_param_dict):
        test_dict["query_attributes"] = query_param_dict["query_attributes"];
    if('segment_size' in query_param_dict):
        test_dict['segment_size'] = query_param_dict['segment_size'];
    else:
        test_dict['segment_size'] = default_segment_size;
    if('produce_GT_field' in query_param_dict):
        test_dict['produce_GT_field'] = query_param_dict['produce_GT_field'];
    if('produce_FILTER_field' in query_param_dict):
        test_dict['produce_FILTER_field'] = query_param_dict['produce_FILTER_field'];
    if('query_block_size' in query_param_dict):
        test_dict['query_block_size'] = query_param_dict['query_block_size'];
    if('query_block_size_margin' in query_param_dict):
        test_dict['query_block_size_margin'] = query_param_dict['query_block_size_margin'];
    if('vid_mapping_file' in test_dict):
        test_dict['vid_mapping_file'] = test_dir+os.path.sep+test_dict['vid_mapping_file'];
    if('callset_mapping_file' in test_dict):
        test_dict['callset_mapping_file'] = test_dir+os.path.sep+test_dict['callset_mapping_file'];
    if('vcf_header_filename' in test_dict):
        for i,val in enumerate(test_dict['vcf_header_filename']):
            test_dict['vcf_header_filename'][i] = test_dir+os.path.sep+val;
    if('reference_genome' in test_dict):
        test_dict['reference_genome'] = test_dir+os.path.sep+test_dict['reference_genome'];
    return test_dict;


loader_json_template_string="""
{
    "row_based_partitioning" : false,
    "column_partitions" : [
        {"begin": 0, "workspace":"", "array": "" }
    ],
    "callset_mapping_file" : "",
    "vid_mapping_file" : "inputs/vid.json",
    "size_per_column_partition": 700 ,
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
    "compress_tiledb_array" : true,
    "segment_size" : 1048576,
    "num_cells_per_tile" : 3
}""";

def create_loader_json(ws_dir, test_name, test_params_dict, col_part, test_dir):
    test_dict=json.loads(loader_json_template_string);
    test_dict['column_partitions'] = col_part;
    for col_part in test_dict['column_partitions']:
        col_part["workspace"] = ws_dir;
        col_part["array"] = test_name+col_part["array"];
    test_dict["callset_mapping_file"] = test_params_dict['callset_mapping_file'];
    if('vid_mapping_file' in test_params_dict):
        test_dict['vid_mapping_file'] = test_params_dict['vid_mapping_file'];
    if('size_per_column_partition' in test_params_dict):
        test_dict['size_per_column_partition'] = test_params_dict['size_per_column_partition'];
    if('segment_size' in test_params_dict):
        test_dict['segment_size'] = test_params_dict['segment_size'];
    else:
        test_dict['segment_size'] = default_segment_size;
    test_dict['vid_mapping_file'] = test_dir+os.path.sep+test_dict['vid_mapping_file'];
    test_dict['callset_mapping_file'] = test_dir+os.path.sep+test_dict['callset_mapping_file'];
    test_dict['vcf_header_filename'] = test_dir+os.path.sep+test_dict['vcf_header_filename'];
    test_dict['reference_genome'] = test_dir+os.path.sep+test_dict['reference_genome'];
    return test_dict;

def add_hdfs_to_loader_json(test_dict, namenode):
    for col_part in test_dict['column_partitions']:
        col_part['workspace'] = namenode+col_part['workspace'];
    return test_dict;

def move_arrays_to_hdfs(ws_dir, namenode):
#    pid = subprocess.Popen('hadoop fs -rm -r '+namenode+ws_dir+'/*', shell=True, stdout=subprocess.PIPE);
#    stdout_string = pid.communicate()[0]
#    if(pid.returncode != 0):
#        sys.stderr.write('Error deleting arrays from workspace in HDFS:'namenode+ws_dir+'\n');
#        sys.exit(-1);
    pid = subprocess.Popen('hadoop fs -put '+ws_dir+'/* '+namenode+ws_dir, shell=True, stdout=subprocess.PIPE);
    stdout_string = pid.communicate()[0]
    if(pid.returncode != 0):
        sys.stderr.write('Error copying array to HDFS workspace:'+namenode+ws_dir+'\n');
        sys.exit(-1);

def get_file_content_and_md5sum(filename):
    with open(filename, 'rb') as fptr:
        data = fptr.read();
        data_list = data.splitlines(True);
        data_list_filter = [k for k in data_list if not k.startswith('##')];
        data_filter = "".join(data_list_filter);
        md5sum_hash_str = str(hashlib.md5(data_filter).hexdigest());

        fptr.close();
        return (data_filter, md5sum_hash_str);

def print_diff(golden_output, test_output):
    print("=======Golden output:=======");
    print(golden_output);
    print("=======Test output:=======");
    print(test_output);
    print("=======END=======");

def cleanup_and_exit(namenode, tmpdir, exit_code):
    if(exit_code == 0):
        shutil.rmtree(tmpdir, ignore_errors=True)
        if("://" in namenode):
            pid = subprocess.Popen('hadoop fs -rm -r '+namenode+tmpdir, shell=True, stdout=subprocess.PIPE);
            pid = subprocess.Popen('hadoop fs -rm -r '+namenode+'/home/hadoop/.tiledb/', shell=True, stdout=subprocess.PIPE);
    sys.exit(exit_code);

def main():
    if(len(sys.argv) < 8):
        sys.stderr.write('Needs 7 arguments <build_dir> <install_dir> <spark_master> <hdfs_namenode> <spark_deploy> <genomicsdb_version> <test_dir>\n');
        sys.exit(-1);
    exe_path = sys.argv[2]+os.path.sep+'bin';
    spark_master = sys.argv[3];
    namenode = sys.argv[4];
    jar_dir = sys.argv[1]+os.path.sep+'target';
    spark_deploy = sys.argv[5];
    genomicsdb_version = sys.argv[6];
    test_dir = sys.argv[7];
    #Switch to tests directory
    parent_dir=os.path.dirname(os.path.realpath(__file__))
    os.chdir(parent_dir)
    hostfile_path=parent_dir+os.path.sep+'hostfile';
    template_vcf_header_path=parent_dir+os.path.sep+'inputs'+os.path.sep+'template_vcf_header.vcf';
    tmpdir = tempfile.mkdtemp()
    ws_dir=tmpdir+os.path.sep+'ws';
    loader_tests = [
            { "name" : "t0_1_2", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                "column_partitions": [
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test0"} ],
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test1"},
                      {"begin": 10000, "workspace":"/tmp/ws", "array": "test2"}
                       ],
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test3"},
                      {"begin": 3000, "workspace":"/tmp/ws", "array": "test4"},
                      {"begin": 6000, "workspace":"/tmp/ws", "array": "test5"},
                      {"begin": 9000, "workspace":"/tmp/ws", "array": "test6"},
                      {"begin": 12000, "workspace":"/tmp/ws", "array": "test7"}
                       ]
                    ],
                "query_params": [
                    { "query_column_ranges" : [12100, 12200], "golden_output": {
                        "spark"   : "golden_outputs/spark_t0_1_2_vcf_at_12100",
                        } },
                    { "query_column_ranges" : [0, 100000], "golden_output": {
                        "spark"   : "golden_outputs/spark_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges" : [12150, 100000], "golden_output": {
                        "spark"   : "golden_outputs/spark_t0_1_2_vcf_at_12150",
                        } },
                    ]
            },
            { "name" : "t0_overlapping", 'golden_output': 'golden_outputs/t0_overlapping',
                'callset_mapping_file': 'inputs/callsets/t0_overlapping.json',
                "column_partitions": [
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test0"} ],
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test1"},
                      {"begin": 10000, "workspace":"/tmp/ws", "array": "test2"}
                       ],
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test3"},
                      {"begin": 3000, "workspace":"/tmp/ws", "array": "test4"},
                      {"begin": 6000, "workspace":"/tmp/ws", "array": "test5"},
                      {"begin": 9000, "workspace":"/tmp/ws", "array": "test6"},
                      {"begin": 12000, "workspace":"/tmp/ws", "array": "test7"}
                       ]
                    ],
                "query_params": [
                    { "query_column_ranges" : [12202, 100000], "golden_output": {
                        "spark"        : "golden_outputs/spark_t0_overlapping_at_12202",
                        }
                    }
                ]
            },
            { "name" : "t6_7_8", 'golden_output' : 'golden_outputs/t6_7_8_loading',
                'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                "column_partitions": [
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test0"} ],
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test1"},
                      {"begin": 500000, "workspace":"/tmp/ws", "array": "test2"},
                      {"begin": 1000000, "workspace":"/tmp/ws", "array": "test3"}
                       ],
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test4"},
                      {"begin": 250000, "workspace":"/tmp/ws", "array": "test5"},
                      {"begin": 500000, "workspace":"/tmp/ws", "array": "test6"},
                      {"begin": 750000, "workspace":"/tmp/ws", "array": "test7"},
                      {"begin": 1000000, "workspace":"/tmp/ws", "array": "test8"}
                       ]
                    ],
                "query_params": [
                    { "query_column_ranges" : [0, 10000000], "golden_output": {
                        "spark": "golden_outputs/spark_t6_7_8_vcf_at_0",
                        },
                      "query_block_size" : 1000000, "query_block_size_margin": 50000 },
                    { "query_column_ranges" : [8029500, 10000000], "golden_output": {
                        "spark": "golden_outputs/spark_t6_7_8_vcf_at_8029500",
                        },
                      "query_block_size" : 100000, "query_block_size_margin": 5000 },
                    { "query_column_ranges" : [8029500, 8029500], "golden_output": {
                        "spark"        : "golden_outputs/spark_t6_7_8_vcf_at_8029500-8029500",
                        } }
                    ]
            },
            { "name" : "t0_1_2_combined", 'golden_output' : 'golden_outputs/t0_1_2_combined',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_combined.json',
                "column_partitions": [
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test0"} ],
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test1"},
                      {"begin": 10000, "workspace":"/tmp/ws", "array": "test2"}
                       ],
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test3"},
                      {"begin": 3000, "workspace":"/tmp/ws", "array": "test4"},
                      {"begin": 6000, "workspace":"/tmp/ws", "array": "test5"},
                      {"begin": 9000, "workspace":"/tmp/ws", "array": "test6"},
                      {"begin": 12000, "workspace":"/tmp/ws", "array": "test7"}
                       ]
                    ],
                "query_params": [
                    { "query_column_ranges" : [0, 1000000], "golden_output": {
                        "spark": "golden_outputs/spark_t0_1_2_combined",
                        },
                      "query_block_size" : 100000, "query_block_size_margin": 5000 },
                    ]
            }, 
            { "name" : "t0_haploid_triploid_1_2_3_triploid_deletion",
                'golden_output' : 'golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_loading',
                'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                'size_per_column_partition': 1200,
                'segment_size': 100,
                "column_partitions": [
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test0"} ],
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test1"},
                      {"begin": 10000, "workspace":"/tmp/ws", "array": "test2"}
                       ],
                    [ {"begin": 0, "workspace":"/tmp/ws", "array": "test3"},
                      {"begin": 3000, "workspace":"/tmp/ws", "array": "test4"},
                      {"begin": 6000, "workspace":"/tmp/ws", "array": "test5"},
                      {"begin": 9000, "workspace":"/tmp/ws", "array": "test6"},
                      {"begin": 12000, "workspace":"/tmp/ws", "array": "test7"}
                       ]
                    ],
                "query_params": [
                    { "query_column_ranges" : [0, 1000000],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'segment_size': 100,
                        "golden_output": {
                        "spark"   : "golden_outputs/spark_t0_haploid_triploid_1_2_3_triploid_deletion_java_vcf",
                        },
                      "query_block_size" : 100000, "query_block_size_margin": 5000 },
                    { "query_column_ranges" : [0, 1000000],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'produce_GT_field': True,
                      'segment_size': 100,
                        "golden_output": {
                        "spark"   : "golden_outputs/spark_t0_haploid_triploid_1_2_3_triploid_deletion_java_vcf_produce_GT",
                        },
                      "query_block_size" : 100000, "query_block_size_margin": 5000 }
                ]
            },
    ];
    if("://" in namenode):
        pid = subprocess.Popen('hadoop fs -mkdir -p '+namenode+'/home/hadoop/.tiledb/', shell=True, stdout=subprocess.PIPE);
        stdout_string = pid.communicate()[0]
        if(pid.returncode != 0):
            sys.stderr.write('Error creating hdfs:///home/hadoop/.tiledb/');
            sys.exit(-1);
    for test_params_dict in loader_tests:
        test_name = test_params_dict['name']
        
        for col_part in test_params_dict['column_partitions']:
            test_loader_dict = create_loader_json(ws_dir, test_name, test_params_dict, col_part, test_dir);
            if(test_name == "t0_1_2"):
                test_loader_dict["compress_tiledb_array"] = True;
            if("://" in namenode):
                test_loader_dict = add_hdfs_to_loader_json(test_loader_dict, namenode);
            loader_json_filename = tmpdir+os.path.sep+test_name+'-loader.json'
            with open(loader_json_filename, 'wb') as fptr:
                json.dump(test_loader_dict, fptr, indent=4, separators=(',', ': '));
                fptr.close();
            # invoke vcf2tiledb -r <rank> where <rank> goes from 0 to num partitions
            # otherwise this only loads the first partition
            for i in range(0, len(col_part)):
                etl_cmd=exe_path+os.path.sep+'vcf2tiledb -r '+str(i)+' '+loader_json_filename
                pid = subprocess.Popen(etl_cmd, shell=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE);
                stdout_string, stderr_string = pid.communicate()
                if(pid.returncode != 0):
                    sys.stderr.write('Loading failed for test: '+test_name+' rank '+str(i)+'\n');
                    sys.stderr.write('Loading command: '+etl_cmd+'\n');
                    sys.stderr.write('Loader file :'+str(test_loader_dict)+'\n');
                    sys.stderr.write('Loading stdout: '+stdout_string+'\n');
                    sys.stderr.write('Loading stderr: '+stderr_string+'\n');
                    cleanup_and_exit(namenode, tmpdir, -1);
            with open(loader_json_filename, 'wb') as fptr:
                json.dump(test_loader_dict, fptr, indent=4, separators=(',', ': '));
                fptr.close();
            for query_param_dict in test_params_dict['query_params']:
                if("://" in namenode):
                    test_query_dict = create_query_json(namenode+ws_dir, test_name, query_param_dict, test_dir)
                else:
                    test_query_dict = create_query_json(ws_dir, test_name, query_param_dict, test_dir)
                test_query_dict['query_attributes'] = vcf_query_attributes_order;
                query_json_filename = tmpdir+os.path.sep+test_name+'-query.json'
                with open(query_json_filename, 'wb') as fptr:
                    json.dump(test_query_dict, fptr, indent=4, separators=(',', ': '));
                    fptr.close();
                loader_argument = loader_json_filename;
                spark_cmd = 'spark-submit --class TestGenomicsDBSparkHDFS --master '+spark_master+' --deploy-mode '+spark_deploy+' --total-executor-cores 1 --executor-memory 512M --conf "spark.yarn.executor.memoryOverhead=3700" --jars '+jar_dir+'/genomicsdb-'+genomicsdb_version+'-jar-with-dependencies.jar '+jar_dir+'/genomicsdb-'+genomicsdb_version+'-examples.jar --loader '+loader_json_filename+' --query '+query_json_filename+' --hostfile '+hostfile_path+' --template_vcf_header '+template_vcf_header_path+' --spark_master '+spark_master+' --jar_dir '+jar_dir;
                pid = subprocess.Popen(spark_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
                stdout_string, stderr_string = pid.communicate()
                if(pid.returncode != 0):
                    sys.stderr.write('Query test: '+test_name+' with query file '+query_json_filename+' failed\n');
                    sys.stderr.write('Spark command was: '+spark_cmd+'\n');
                    sys.stderr.write('Spark stdout was: '+stdout_string+'\n');
                    sys.stderr.write('Spark stderr was: '+stderr_string+'\n');
                    cleanup_and_exit(namenode, tmpdir, -1);
                stdout_list = stdout_string.splitlines(True);
                stdout_list_filter = [k for k in stdout_list if not k.startswith('##')];
                stdout_filter = "".join(stdout_list_filter);
                md5sum_hash_str = str(hashlib.md5(stdout_filter).hexdigest())
                if('golden_output' in query_param_dict and 'spark' in query_param_dict['golden_output']):
                    golden_stdout, golden_md5sum = get_file_content_and_md5sum(query_param_dict['golden_output']['spark']);
                    if(golden_md5sum != md5sum_hash_str):
                        sys.stdout.write('Mismatch in query test: '+test_name+' with column ranges: '+str(query_param_dict['query_column_ranges'])+' and loaded with '+str(len(col_part))+' partitions\n');
                        print_diff(golden_stdout, stdout_filter);
                        sys.stderr.write('Spark stdout was: '+stdout_string+'\n');
                        sys.stderr.write('Spark stderr was: '+stderr_string+'\n');
                        cleanup_and_exit(namenode, tmpdir, -1);
                    else:
                        sys.stdout.write('Query test: '+test_name+' with column ranges: '+str(query_param_dict['query_column_ranges'])+' and loaded with '+str(len(col_part))+' partitions passed\n');
    cleanup_and_exit(namenode, tmpdir, 0); 

if __name__ == '__main__':
    main()
