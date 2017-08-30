#!/usr/bin/env python

#The MIT License (MIT)
#Copyright (c) 2016-2017 Intel Corporation

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
from collections import OrderedDict
import jsondiff

query_json_template_string="""
{   
        "workspace" : "",
        "array" : "",
        "vcf_header_filename" : ["inputs/template_vcf_header.vcf"],
        "query_column_ranges" : [ [ [0, 10000000000 ] ] ],
        "query_row_ranges" : [ [ [0, 3 ] ] ],
        "reference_genome" : "inputs/chr1_10MB.fasta.gz",
        "query_attributes" : [ "REF", "ALT", "BaseQRankSum", "MQ", "RAW_MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP", "PID", "PGT" ]
}"""

vcf_query_attributes_order = [ "END", "REF", "ALT", "BaseQRankSum", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "MQ", "RAW_MQ", "MQ0", "DP", "GT", "GQ", "SB", "AD", "PL", "PGT", "PID", "MIN_DP", "DP_FORMAT", "FILTER" ];
query_attributes_with_DS_ID = [ "REF", "ALT", "BaseQRankSum", "MQ", "RAW_MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP", "PID", "PGT", "DS", "ID" ];
query_attributes_with_PL_only = [ "PL" ]
query_attributes_with_MLEAC_only = [ "MLEAC" ]
default_segment_size = 40

def create_query_json(ws_dir, test_name, query_param_dict):
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
    "compress_tiledb_array" : false,
    "segment_size" : 1048576,
    "num_cells_per_tile" : 3
}""";

def create_loader_json(ws_dir, test_name, test_params_dict):
    test_dict=json.loads(loader_json_template_string);
    if('column_partitions' in test_params_dict):
        test_dict['column_partitions'] = test_params_dict['column_partitions'];
    test_dict["column_partitions"][0]["workspace"] = ws_dir;
    test_dict["column_partitions"][0]["array"] = test_name;
    test_dict["callset_mapping_file"] = test_params_dict['callset_mapping_file'];
    if('vid_mapping_file' in test_params_dict):
        test_dict['vid_mapping_file'] = test_params_dict['vid_mapping_file'];
    if('size_per_column_partition' in test_params_dict):
        test_dict['size_per_column_partition'] = test_params_dict['size_per_column_partition'];
    if('segment_size' in test_params_dict):
        test_dict['segment_size'] = test_params_dict['segment_size'];
    else:
        test_dict['segment_size'] = default_segment_size;
    return test_dict;

def get_file_content_and_md5sum(filename):
    with open(filename, 'rb') as fptr:
        data = fptr.read();
        md5sum_hash_str = str(hashlib.md5(data).hexdigest())
        fptr.close();
        return (data, md5sum_hash_str);

def print_diff(golden_output, test_output):
    print("=======Golden output:=======");
    print(golden_output);
    print("=======Test output:=======");
    print(test_output);
    print("=======END=======");

def cleanup_and_exit(tmpdir, exit_code):
    if(exit_code == 0):
        shutil.rmtree(tmpdir, ignore_errors=True)
    sys.exit(exit_code);

def main():
    #lcov gcda directory prefix
    gcda_prefix_dir = '../';
    if(len(sys.argv) < 3):
        sys.stderr.write('Needs 2 arguments <build_dir> <install_dir>\n');
        sys.exit(-1);
    gcda_prefix_dir = sys.argv[1];
    exe_path = sys.argv[2]+os.path.sep+'bin';
    #Switch to tests directory
    parent_dir=os.path.dirname(os.path.realpath(__file__))
    os.chdir(parent_dir)
    #Zero line coverage
    subprocess.call('lcov --directory '+gcda_prefix_dir+' --zerocounters', shell=True);
    tmpdir = tempfile.mkdtemp()
    ws_dir=tmpdir+os.path.sep+'ws';
    loader_tests = [
            { "name" : "t0_1_2", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges" : [0, 1000000000],
                        "query_attributes": query_attributes_with_PL_only,
                        "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0_with_PL_only",
                        } },
                    { "query_column_ranges" : [0, 1000000000],#vid and callset jsons passed through query json
                        "query_without_loader": True,
                        "vid_mapping_file": "inputs/vid.json",
                        "callset_mapping_file": "inputs/callsets/t0_1_2.json",
                        "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges" : [12150, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12150",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_12150",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } },
                    { "query_column_ranges" : [0, 1000000000],
                        "produce_FILTER_field": True, "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0_with_FILTER",
                        } },
                    ]
            },
            { "name" : "t0_1_2_csv", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_csv.json',
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges" : [12150, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12150",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_12150",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } }
                    ]
            },
            { "name" : "t0_overlapping", 'golden_output': 'golden_outputs/t0_overlapping',
                'callset_mapping_file': 'inputs/callsets/t0_overlapping.json',
                "query_params": [
                    { "query_column_ranges" : [12202, 1000000000], "golden_output": {
                        "vcf"        : "golden_outputs/t0_overlapping_at_12202",
                        }
                    }
                ]
            },
            { "name" : "t0_overlapping_at_12202", 'golden_output': 'golden_outputs/t0_overlapping_at_12202',
                'callset_mapping_file': 'inputs/callsets/t0_overlapping.json',
                'column_partitions': [ {"begin": 12202, "workspace":"", "array": "" }]
            },
            { "name" : "t6_7_8", 'golden_output' : 'golden_outputs/t6_7_8_loading',
                'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_0",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_0",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_0",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_0",
                        } },
                    { "query_column_ranges" : [8029500, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_8029500",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_8029500",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_8029500",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_8029500",
                        } },
                    { "query_column_ranges" : [8029500, 8029500], "golden_output": {
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_8029500-8029500",
                        } }
                    ]
            },
            { "name" : "java_t0_1_2", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0_phased_GT",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0_phased_GT",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges" : [12150, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12150_phased_GT",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_12150_phased_GT",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } }
                    ]
            },
            { "name" : "java_buffer_stream_t0_1_2", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_buffer.json',
                'stream_name_to_filename_mapping': 'inputs/callsets/t0_1_2_buffer_mapping.json',
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges" : [12150, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12150",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_12150",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } }
                    ]
            },
            { "name" : "java_buffer_stream_multi_contig_t0_1_2", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_buffer.json',
                'stream_name_to_filename_mapping': 'inputs/callsets/t0_1_2_buffer_mapping.json',
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges" : [12150, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12150",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_12150",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } }
                    ]
            },
            { "name" : "test_new_fields", 'golden_output' : 'golden_outputs/t6_7_8_new_field_gatk.vcf',
                'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                'vid_mapping_file': 'inputs/vid_MLEAC_MLEAF.json',
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000],
                      "query_attributes" : query_attributes_with_MLEAC_only, "golden_output": {
                        "calls"        : "golden_outputs/test_new_fields_MLEAC_only.json",
                        } },
                    ]
            },
            { "name" : "test_info_combine_ops0", 'golden_output' : 'golden_outputs/info_ops0.vcf',
                'callset_mapping_file': 'inputs/callsets/info_ops.json',
                'vid_mapping_file': 'inputs/vid_info_ops0.json'
            },
            { "name" : "test_info_combine_ops1", 'golden_output' : 'golden_outputs/info_ops1.vcf',
                'callset_mapping_file': 'inputs/callsets/info_ops.json',
                'vid_mapping_file': 'inputs/vid_info_ops1.json'
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t0_1_2",
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                'chromosome_interval': '1:1-100000000',
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000],
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges" : [12150, 1000000000],
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } }
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t6_7_8",
                'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                'chromosome_interval': '1:1-100000000',
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000],
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                        "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_0_phased_GT",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_0_phased_GT",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_0",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_0",
                        } },
                    { "query_column_ranges" : [8029500, 1000000000],
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                        "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_8029500_phased_GT",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_8029500_phased_GT",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_8029500",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_8029500",
                        } }
                    ]
            },
            { "name" : "t0_1_2_combined", 'golden_output' : 'golden_outputs/t0_1_2_combined',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_combined.json',
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000], "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_combined",
                        "batched_vcf": "golden_outputs/t0_1_2_combined",
                        } },
                    ]
            }, 
            { "name" : "test_flag_field", 'golden_output' : 'golden_outputs/t0_1_2_DS_ID_vcf_at_0',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                'vid_mapping_file': 'inputs/vid_DS_ID.json',
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000],
                        "query_attributes": query_attributes_with_DS_ID, "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_DS_ID_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_DS_ID_variants_at_0",
                        } },
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t0_1_2_with_DS_ID",
                'vid_mapping_file': 'inputs/vid_DS_ID_phased_GT.json',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                'chromosome_interval': '1:1-100000000',
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000],
                        'vid_mapping_file': 'inputs/vid_DS_ID_phased_GT.json',
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "query_attributes": query_attributes_with_DS_ID, "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_DS_ID_calls_at_0_phased_GT",
                        "variants"   : "golden_outputs/t0_1_2_DS_ID_variants_at_0_phased_GT",
                        } },
                    ]
            },
            { "name" : "t0_1_2_as_array", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_as_array.json',
                "vid_mapping_file": "inputs/vid_as_array.json",
            },
            { "name" : "t0_with_missing_PL_SB_fields", 'golden_output' : 'golden_outputs/t0_with_missing_PL_SB_fields_t1.vcf',
                'callset_mapping_file': 'inputs/callsets/t0_with_missing_PL_SB_fields_t1.json',
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000], "golden_output": {
                        "calls"      : "golden_outputs/t0_with_missing_PL_SB_fields_t1_calls.json",
                        } },
                    ]
            },
            { "name" : "t0_haploid_triploid_1_2_3_triploid_deletion",
                'golden_output' : 'golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_loading',
                'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                'size_per_column_partition': 1200,
                'segment_size': 100,
                "query_params": [
                    { "query_column_ranges" : [0, 1000000000],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'segment_size': 100,
                        "golden_output": {
                        "vcf"        : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_vcf",
                        "java_vcf"   : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_java_vcf",
                        } },
                    { "query_column_ranges" : [0, 1000000000],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'produce_GT_field': True,
                      'segment_size': 100,
                        "golden_output": {
                        "vcf"        : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_vcf_produce_GT",
                        "java_vcf"   : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_java_vcf_produce_GT",
                        } }
                ]
            },
    ];
    for test_params_dict in loader_tests:
        test_name = test_params_dict['name']
        test_loader_dict = create_loader_json(ws_dir, test_name, test_params_dict);
        if(test_name == "t0_1_2"):
            test_loader_dict["compress_tiledb_array"] = True;
        loader_json_filename = tmpdir+os.path.sep+test_name+'.json'
        with open(loader_json_filename, 'wb') as fptr:
            json.dump(test_loader_dict, fptr, indent=4, separators=(',', ': '));
            fptr.close();
        if(test_name  == 'java_t0_1_2'):
            pid = subprocess.Popen('java -ea TestGenomicsDB --load '+loader_json_filename, shell=True,
                    stdout=subprocess.PIPE);
        elif(test_name == 'java_buffer_stream_multi_contig_t0_1_2'):
            pid = subprocess.Popen('java -ea TestBufferStreamGenomicsDBImporter -iterators '+loader_json_filename+' '
                    +test_params_dict['stream_name_to_filename_mapping']
                    +' 1024 0 0 100 true ',
                    shell=True, stdout=subprocess.PIPE);
        elif(test_name == 'java_buffer_stream_t0_1_2'):
            pid = subprocess.Popen('java -ea TestBufferStreamGenomicsDBImporter '+loader_json_filename
                    +' '+test_params_dict['stream_name_to_filename_mapping'],
                    shell=True, stdout=subprocess.PIPE);
        elif(test_name.find('java_genomicsdb_importer_from_vcfs') != -1):
            arg_list = ' -L '+test_params_dict['chromosome_interval'] + ' -w ' + ws_dir + ' -A '+test_name \
                    +' --use_samples_in_order ' + ' --batchsize=2 ';
            with open(test_params_dict['callset_mapping_file'], 'rb') as cs_fptr:
                callset_mapping_dict = json.load(cs_fptr, object_pairs_hook=OrderedDict)
                for callset_name, callset_info in callset_mapping_dict['callsets'].iteritems():
                    arg_list += ' '+callset_info['filename'];
                cs_fptr.close();
            pid = subprocess.Popen('java -ea TestGenomicsDBImporterWithMergedVCFHeader '+arg_list,
                    shell=True, stdout=subprocess.PIPE);
        else:
            pid = subprocess.Popen(exe_path+os.path.sep+'vcf2tiledb '+loader_json_filename, shell=True,
                    stdout=subprocess.PIPE);
        stdout_string = pid.communicate()[0]
        if(pid.returncode != 0):
            sys.stderr.write('Loader test: '+test_name+' failed\n');
            cleanup_and_exit(tmpdir, -1);
        md5sum_hash_str = str(hashlib.md5(stdout_string).hexdigest())
        if('golden_output' in test_params_dict):
            golden_stdout, golden_md5sum = get_file_content_and_md5sum(test_params_dict['golden_output']);
            if(golden_md5sum != md5sum_hash_str):
                sys.stderr.write('Loader stdout mismatch for test: '+test_name+'\n');
                print_diff(golden_stdout, stdout_string);
                cleanup_and_exit(tmpdir, -1);
        if('query_params' in test_params_dict):
            for query_param_dict in test_params_dict['query_params']:
                test_query_dict = create_query_json(ws_dir, test_name, query_param_dict)
                query_types_list = [
                        ('calls','--print-calls'),
                        ('variants',''),
                        ('vcf','--produce-Broad-GVCF'),
                        ('batched_vcf','--produce-Broad-GVCF -p 128'),
                        ('java_vcf', ''),
                        ('consolidate_and_vcf', '--produce-Broad-GVCF'), #keep as the last query test
                        ]
                for query_type,cmd_line_param in query_types_list:
                    if(query_type == 'vcf' or query_type == 'batched_vcf' or query_type.find('java_vcf') != -1):
                        test_query_dict['query_attributes'] = vcf_query_attributes_order;
                    query_json_filename = tmpdir+os.path.sep+test_name+'_'+query_type+'.json'
                    with open(query_json_filename, 'wb') as fptr:
                        json.dump(test_query_dict, fptr, indent=4, separators=(',', ': '));
                        fptr.close();
                    query_command = ''
                    if(query_type == 'java_vcf'):
                        loader_argument = loader_json_filename;
                        if("query_without_loader" in query_param_dict and query_param_dict["query_without_loader"]):
                            loader_argument = '""'
                        query_command = 'java -ea TestGenomicsDB --query -l '+loader_argument+' '+query_json_filename;
                        pid = subprocess.Popen(query_command, shell=True, stdout=subprocess.PIPE);
                    else:
                        if(query_type == 'consolidate_and_vcf'):
                            retcode = subprocess.call(exe_path+os.path.sep+'consolidate_tiledb_array '+ws_dir+' '+test_name,
                                    shell=True)
                            if(retcode != 0):
                                sys.stderr.write('TileDB array consolidation failed '+ws_dir+' '+test_name+'\n');
                                cleanup_and_exit(tmpdir, -1);
                        loader_argument = ' -l '+loader_json_filename;
                        if("query_without_loader" in query_param_dict and query_param_dict["query_without_loader"]):
                            loader_argument = ''
                        query_command = (exe_path+os.path.sep+'gt_mpi_gather -s %d'+loader_argument
                            + ' -j '
                            +query_json_filename+' '+cmd_line_param)%(test_query_dict['segment_size']);
                        pid = subprocess.Popen(query_command, shell=True, stdout=subprocess.PIPE);
                    stdout_string = pid.communicate()[0]
                    if(pid.returncode != 0):
                        sys.stderr.write('Command '+query_command+'\n')
                        sys.stderr.write('Query test: '+test_name+'-'+query_type+' failed\n');
                        cleanup_and_exit(tmpdir, -1);
                    md5sum_hash_str = str(hashlib.md5(stdout_string).hexdigest())
                    if('golden_output' in query_param_dict and query_type in query_param_dict['golden_output']):
                        golden_stdout, golden_md5sum = get_file_content_and_md5sum(query_param_dict['golden_output'][query_type]);
                        if(golden_md5sum != md5sum_hash_str):
                            is_error = True;
                            #do JSON diff for variant and call format print
                            json_diff_result = None
                            if(query_type in set(['calls', 'variants'])):
                                try:
                                    golden_stdout_dict = json.loads(golden_stdout);
                                    test_stdout_dict = json.loads(stdout_string);
                                    json_diff_result = jsondiff.diff(golden_stdout_dict, test_stdout_dict);
                                    if(len(json_diff_result) == 0):
                                        is_error = False;
                                except:
                                    json_diff_result = None;
                                    is_error = True
                            if(is_error):
                                sys.stderr.write('Mismatch in query test: '+test_name+'-'+query_type+'\n');
                                print_diff(golden_stdout, stdout_string);
                                if(json_diff_result):
                                    print(json.dumps(json_diff_result, indent=4, separators=(',', ': ')));
                                cleanup_and_exit(tmpdir, -1);
    coverage_file='coverage.info'
    subprocess.call('lcov --directory '+gcda_prefix_dir+' --capture --output-file '+coverage_file, shell=True);
    #Remove protocol buffer generated files from the coverage information
    subprocess.call("lcov --remove "+coverage_file+" '/opt*' '/usr*' 'dependencies*' '*.pb.h' '*.pb.cc' -o "+coverage_file, shell=True);
    cleanup_and_exit(tmpdir, 0); 

if __name__ == '__main__':
    main()
