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
        "array_name" : "",
        "vcf_header_filename" : ["inputs/template_vcf_header.vcf"],
        "query_column_ranges": [{
            "range_list": [{
                "low": 0,
                "high": 10000000000
            }]
        }],
        "query_row_ranges": [{
            "range_list": [{
                "low": 0,
                "high": 3
            }]
        }],
        "reference_genome" : "inputs/chr1_10MB.fasta.gz",
        "attributes" : [ "REF", "ALT", "BaseQRankSum", "MQ", "RAW_MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP", "PID", "PGT" ]
}"""

vcf_attributes_order = [ "END", "REF", "ALT", "BaseQRankSum", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "MQ", "RAW_MQ", "MQ0", "DP", "GT", "GQ", "SB", "AD", "PL", "PGT", "PID", "MIN_DP", "DP_FORMAT", "FILTER" ];
asa_vcf_attributes = [ "END", "REF", "ALT", "BaseQRankSum", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "MQ", "RAW_MQ", "MQ0", "DP", "GT", "GQ", "SB", "AD", "PL", "PGT", "PID", "MIN_DP", "DP_FORMAT", "FILTER", "AS_RAW_MQ", "AS_RAW_MQRankSum" ];
attributes_with_DS_ID = [ "REF", "ALT", "BaseQRankSum", "MQ", "RAW_MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP", "PID", "PGT", "DS", "ID" ];
attributes_with_PL_only = [ "PL" ]
attributes_with_MLEAC_only = [ "MLEAC" ]
default_segment_size = 40

def create_query_json(ws_dir, test_name, query_param_dict):
    test_dict=json.loads(query_json_template_string);
    test_dict["workspace"] = ws_dir
    test_dict["array_name"] = test_name
    if('query_column_ranges' in query_param_dict):
        test_dict["query_column_ranges"] = query_param_dict["query_column_ranges"]
    else:
        test_dict['scan_full'] = True
    if("vid_mapping_file" in query_param_dict):
        test_dict["vid_mapping_file"] = query_param_dict["vid_mapping_file"];
    if("callset_mapping_file" in query_param_dict):
        test_dict["callset_mapping_file"] = query_param_dict["callset_mapping_file"];
    if("attributes" in query_param_dict):
        test_dict["attributes"] = query_param_dict["attributes"];
    if('segment_size' in query_param_dict):
        test_dict['segment_size'] = query_param_dict['segment_size'];
    else:
        test_dict['segment_size'] = default_segment_size;
    if('produce_GT_field' in query_param_dict):
        test_dict['produce_GT_field'] = query_param_dict['produce_GT_field'];
    if('produce_FILTER_field' in query_param_dict):
        test_dict['produce_FILTER_field'] = query_param_dict['produce_FILTER_field'];
    if('sites_only_query' in query_param_dict):
        test_dict['sites_only_query'] = query_param_dict['sites_only_query']
    if('produce_GT_with_min_PL_value_for_spanning_deletions' in query_param_dict):
        test_dict['produce_GT_with_min_PL_value_for_spanning_deletions'] = \
                query_param_dict['produce_GT_with_min_PL_value_for_spanning_deletions']
    return test_dict;


loader_json_template_string="""
{
    "row_based_partitioning" : false,
    "column_partitions" : [
        {"begin": 0, "workspace":"", "array_name": "" }
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
    test_dict["column_partitions"][0]["array_name"] = test_name;
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

def modify_query_column_ranges_for_PB(test_query_dict):
    if('query_column_ranges' in test_query_dict):
        original_query_column_ranges = test_query_dict['query_column_ranges']
        new_query_column_ranges = []
        for curr_entry in original_query_column_ranges:
            if(type(curr_entry) is dict and 'range_list' in curr_entry):
                new_interval_list = []
                for curr_interval in curr_entry['range_list']:
                    if(type(curr_interval) is dict and 'low' in curr_interval
                            and 'high' in curr_interval):
                        new_interval_list.append({'column_interval': { 'column_interval':
                            { 'begin': curr_interval['low'], 'end': curr_interval['high'] } } })
                new_entry = { 'column_or_interval_list': new_interval_list }
                new_query_column_ranges.append(new_entry)
        test_query_dict['query_column_ranges'] = new_query_column_ranges



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
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [ [ [0, 1000000] ] ],
                      "pass_through_query_json": True,
                      "golden_output": {
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges" : [{
                        "range_list": [{
                             "low": 0,
                             "high": 1000000000
                          }]
                      }],
                      "sites_only_query": True,
                      "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_vcf_sites_only_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_sites_only_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12100,
                            "high": 12100
                        }]
                    }], "golden_output": {   #
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12100",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12100,
                            "high": 12100
                        },{
                            "low": 12141,
                            "high": 12141
                        }]
                    }], "golden_output": {   #
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12100_12141",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12100,
                            "high": 12100
                        },{
                            "low": 12141,
                            "high": 12141
                        },{
                            "low": 12150,
                            "high": 12150
                        }]
                    }], "golden_output": {   #
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12100_12141_12150",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12100,
                            "high": 12100
                        },
                            {
                                "low": 12141,
                                "high": 12150
                            }]
                    }], "golden_output": {   #
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12100_12141_to_12150",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12100,
                            "high": 12100
                        },
                            {
                                "low": 12141,
                                "high": 12150
                            },
                            {
                                "low": 12300,
                                "high": 12300
                            },
                            {
                                "low": 17384,
                                "high": 17384
                            }]
                    }], "golden_output": {   #
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12100_12141_to_12150_12300_17384",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "attributes": attributes_with_PL_only,
                        "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0_with_PL_only",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],#vid and callset jsons passed through query json
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
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12150",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_12150",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "produce_FILTER_field": True, "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0_with_FILTER",
                        } },
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t0_1_2_multi_contig",
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                'chromosome_intervals': [ '1:1-12160', '1:12161-12200', '1:12201-18000' ],
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                'generate_array_name_from_partition_bounds': True,
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 18000
                        }]
                    }],
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_genomicsdb_importer_from_vcfs_t0_1_2_multi_contig_vcf_0_18000",
                        } },
                    {
                        "query_contig_interval": {
                            "contig": "1",
                            "begin": 12151,
                            "end": 18000
                        },
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_genomicsdb_importer_from_vcfs_t0_1_2_multi_contig_vcf_12150_18000",
                        } }
                ]
            },
            { "name" : "t0_1_2_csv", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_csv.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
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
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12202,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "vcf"        : "golden_outputs/t0_overlapping_at_12202",
                        }
                    }
                ]
            },
            { "name" : "t0_overlapping_at_12202", 'golden_output': 'golden_outputs/t0_overlapping_at_12202',
                'callset_mapping_file': 'inputs/callsets/t0_overlapping.json',
                'column_partitions': [ {"begin": 12202, "workspace":"", "array_name": "" }]
            },
            { "name" : "t6_7_8", 'golden_output' : 'golden_outputs/t6_7_8_loading',
                'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_0",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_0",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_0",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_0",
                        } },
		    { "query_column_ranges" : [{
                        "range_list": [{
                          "low": 0,
                          "high": 1000000000
                        }]
                       }],
		      "sites_only_query": True, "golden_output": {
                        "vcf"        : "golden_outputs/t6_7_8_vcf_sites_only_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 8029500,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_8029500",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_8029500",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_8029500",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_8029500",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 8029500,
                            "high": 8029500
                        }]
                    }], "golden_output": {
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_8029500-8029500",
                        } }
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t6_7_8_multi_contig",
                'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                'chromosome_intervals': [ '1:1-8029500','1:8029501-8029501', '1:8029502-10000000' ],
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                'generate_array_name_from_partition_bounds': True,
                "query_params": [
                    {   'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_t6_7_8_vcf_at_0",
                        } },
                    {
                        "query_contig_interval": {
                            "contig": "1",
                            "begin": 8029501,
                            "end": 8029510
                        },
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_t6_7_8_vcf_at_8029500",
                        } },
                    {
                        "query_contig_interval": {
                            "contig": "1",
                            "begin": 8029502,
                            "end": 8029502
                        },
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_t6_7_8_vcf_at_8029501",
                        } }
                ]
            },
            { "name" : "java_t0_1_2", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0_phased_GT",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0_phased_GT",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
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
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
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
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
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
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                      "attributes" : attributes_with_MLEAC_only, "golden_output": {
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
                'chromosome_intervals': [ '1:1-100000000' ],
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }],
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
                'chromosome_intervals': [ '1:1-100000000' ],
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                        "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_0_phased_GT",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_0_phased_GT",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_0",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 8029500,
                            "high": 1000000000
                        }]
                    }],
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
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_combined",
                        "batched_vcf": "golden_outputs/t0_1_2_combined",
                        } },
                    ]
            }, 
            { "name" : "test_flag_field", 'golden_output' : 'golden_outputs/t0_1_2_DS_ID_vcf_at_0',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                'vid_mapping_file': 'inputs/vid_DS_ID.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "attributes": attributes_with_DS_ID, "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_DS_ID_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_DS_ID_variants_at_0",
                        } },
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t0_1_2_with_DS_ID",
                'vid_mapping_file': 'inputs/vid_DS_ID_phased_GT.json',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                'chromosome_intervals': [ '1:1-100000000' ],
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        'vid_mapping_file': 'inputs/vid_DS_ID_phased_GT.json',
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "attributes": attributes_with_DS_ID, "golden_output": {
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
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
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
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'segment_size': 100,
                      "golden_output": {
                        "vcf"        : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_vcf",
                        "java_vcf"   : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_java_vcf",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                        }],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'produce_GT_field': True,
                      'segment_size': 100,
                      "golden_output": {
                        "vcf"        : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_vcf_produce_GT",
                        "java_vcf"   : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_java_vcf_produce_GT",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                        }],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'produce_GT_field': True,
                      'produce_GT_with_min_PL_value_for_spanning_deletions': True,
                      'segment_size': 100,
                      "golden_output": {
                        "vcf"        : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_vcf_produce_GT_for_min_value_PL",
                        "java_vcf"   : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_java_vcf_produce_GT_for_min_PL",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                        }],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'sites_only_query': True,
                      'segment_size': 100,
                      "golden_output": {
                        "vcf"        : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_vcf_sites_only",
                        "java_vcf"   : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_java_vcf_sites_only",
                        } },
                ]
            },
            { "name" : "t0_1_2_all_asa", 'golden_output' : 'golden_outputs/t0_1_2_all_asa_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_all_asa.json',
                'vid_mapping_file': 'inputs/vid_all_asa.json',
                'size_per_column_partition': 3000,
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                      "force_override": True,
                      'segment_size': 100,
                      "attributes": asa_vcf_attributes,
                        "golden_output": {
                        "vcf"      : "golden_outputs/t0_1_2_all_asa_loading",
                        } },
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t0_1_2_all_asa",
                'callset_mapping_file': 'inputs/callsets/t0_1_2_all_asa.json',
                'vid_mapping_file': 'inputs/vid_all_asa.json',
                'chromosome_intervals': [ '1:1-100000000' ],
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "force_override": True,
                        'segment_size': 100,
                        "attributes": asa_vcf_attributes,
                        "golden_output": {
                        "vcf"      : "golden_outputs/t0_1_2_all_asa_loading",
                        "java_vcf"   : "golden_outputs/t0_1_2_all_asa_java_query_vcf",
                        } },
                    ]
            },
            { "name" : "min_PL_spanning_deletion", 'golden_output' : 'golden_outputs/min_PL_spanning_deletion_load_stdout',
                'callset_mapping_file': 'inputs/callsets/min_PL_spanning_deletion.json',
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                      'produce_GT_field': True, "golden_output": {
                        "vcf"        : "golden_outputs/min_PL_spanning_deletion_vcf_no_min_PL",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                      'produce_GT_field': True,
                      'produce_GT_with_min_PL_value_for_spanning_deletions': True,
                      "golden_output": {
                        "vcf"        : "golden_outputs/min_PL_spanning_deletion_vcf",
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
            arg_list = ''
            for interval in test_params_dict['chromosome_intervals']:
                arg_list += ' -L '+interval
            arg_list += ' -w ' + ws_dir +' --use_samples_in_order ' + ' --batchsize=2 '
            arg_list += ' --vidmap-output '+ tmpdir + os.path.sep + 'vid.json'
            arg_list += ' --callset-output '+ tmpdir + os.path.sep + 'callsets.json'
            if('generate_array_name_from_partition_bounds' not in test_params_dict or
                    not test_params_dict['generate_array_name_from_partition_bounds']):
                arg_list += ' -A ' + test_name
            with open(test_params_dict['callset_mapping_file'], 'rb') as cs_fptr:
                callset_mapping_dict = json.load(cs_fptr, object_pairs_hook=OrderedDict)
                for callset_name, callset_info in callset_mapping_dict['callsets'].iteritems():
                    arg_list += ' '+callset_info['filename'];
                cs_fptr.close();
            pid = subprocess.Popen('java -ea TestGenomicsDBImporterWithMergedVCFHeader --size_per_column_partition 16384 '
                                   '--segment_size 10485760'+arg_list,
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
                if(test_name.find('java_genomicsdb_importer_from_vcfs') != -1 and
                        'generate_array_name_from_partition_bounds' in test_params_dict
                        and test_params_dict['generate_array_name_from_partition_bounds']):
                    if('array' in test_query_dict):
                        del test_query_dict['array']
                    if('array_name' in test_query_dict):
                        del test_query_dict['array_name']
                    if('query_column_ranges' in test_query_dict):
                        del test_query_dict['query_column_ranges']
                        test_query_dict['scan_full'] = True
                query_types_list = [
                        ('calls','--print-calls'),
                        ('variants',''),
                        ('vcf','--produce-Broad-GVCF'),
                        ('batched_vcf','--produce-Broad-GVCF -p 128'),
                        ('java_vcf', ''),
                        ('consolidate_and_vcf', '--produce-Broad-GVCF'), #keep as the last query test
                        ]
                for query_type,cmd_line_param in query_types_list:
                    if('golden_output' in query_param_dict and query_type in query_param_dict['golden_output']):
                        if((query_type == 'vcf' or query_type == 'batched_vcf' or query_type.find('java_vcf') != -1)
                                and 'force_override' not in query_param_dict):
                            test_query_dict['attributes'] = vcf_attributes_order;
                        if(query_type.find('java_vcf') != -1 and 'pass_through_query_json' not in query_param_dict):
                            modify_query_column_ranges_for_PB(test_query_dict)
                        query_json_filename = tmpdir+os.path.sep+test_name+'_'+query_type+'.json'
                        with open(query_json_filename, 'wb') as fptr:
                            json.dump(test_query_dict, fptr, indent=4, separators=(',', ': '));
                            fptr.close();
                        if(query_type == 'java_vcf'):
                            loader_argument = loader_json_filename;
                            misc_args = ''
                            if("query_without_loader" in query_param_dict and query_param_dict["query_without_loader"]):
                                loader_argument = '""'
                            if("pass_through_query_json" in query_param_dict and query_param_dict["pass_through_query_json"]):
                                misc_args = "--pass_through_query_json"
                            if('query_contig_interval' in query_param_dict):
                                query_contig_interval_dict = query_param_dict['query_contig_interval']
                                misc_args += ('--chromosome '+query_contig_interval_dict['contig'] \
                                        + ' --begin %d --end %d')%(query_contig_interval_dict['begin'],
                                                query_contig_interval_dict['end'])
                            query_command = 'java -ea TestGenomicsDB --query -l '+loader_argument+' '+query_json_filename \
                                + ' ' + misc_args;
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
        shutil.rmtree(ws_dir, ignore_errors=True)
    coverage_file='coverage.info'
    subprocess.call('lcov --directory '+gcda_prefix_dir+' --capture --output-file '+coverage_file, shell=True);
    #Remove protocol buffer generated files from the coverage information
    subprocess.call("lcov --remove "+coverage_file+" '/opt*' '/usr*' 'dependencies*' '*.pb.h' '*.pb.cc' -o "+coverage_file, shell=True);
    cleanup_and_exit(tmpdir, 0); 

if __name__ == '__main__':
    main()
