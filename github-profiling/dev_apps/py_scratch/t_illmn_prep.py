#! /usr/bin/python3
'''
DB-139 this script generates interesting_pos_n for density evealuation
  #parallel is 12
1) create loader_cfg 16 partitions 5000, 10000 and 15000 samples (relevant callset and hist)
2) create 12 queries, one for each loaded partitions
   gt_mpi_run -j <query_i.json> --produce-interesting-positions > interesting_pos
OUTPUT @ /path/to/*.json

RUN loader
1. mkdir -p $dest_path
2. touch $dest_path/__tiledb_workspace.tdb
    nohup /usr/bin/time -v mpirun -np 12 vcf2tiledb /path/to/16-lc_12_illmn_prep5k.json > ~/path/to/log 2>&1
'''

import sys
import os
import os.path
from collections import namedtuple
import json

from pygenomicsdblib.utils.histogram import HistogramManager
from pygenomicsdblib.utils.constant_def import QUERY_FULL_SCAN_CFG_T, LOADER_CFG_T, __COL_PARTITION_SIZE_UNIT
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if not root_path in sys.path:
    sys.path.append(root_path)
import shared.profiling_env as ppenv

num_samples = [1000, 3000, 7500]    # [1000, 3000, 5000, 7500, 10000, 15000]

CREATE_COL_PARTITION_F = lambda start_pos, tdbws: [{"array": "TEST%d" % (i), "begin": begin, "workspace": tdbws} \
            for i, begin in enumerate(start_pos)]

num_parallel = 16
num_partition = 16
is_ssd = False

TEST_NAME = "illmn"

# for ref. already have 5K interest_pos
# hist_list = ['histogram.5000', 'histogram.10000', 'hist_illmn_10000.json']

histogram_fpath = os.path.join(ppenv.get_histogram_root(), 'histogram.10000')
array_start_pos = HistogramManager(histogram_fpath).calc_bin_begin_pos(num_partition)
jsonconfig = namedtuple('TestConfig', ["test_name", "callsets", "vid_mapping", \
    "ref_geno", "is_ssd", "array_start_pos", "num_parallel"])

IllmnQueryAttr = ["REF", "ALT", "QUAL", "GT"]

json_cfg = jsonconfig(TEST_NAME, "callsets_illmn_15000.json", "illmn_vid.json", \
    "Homo_sapiens_assembly19.fasta", is_ssd, array_start_pos, num_parallel)

tdb_ws_root = os.path.join(ppenv.get_tdb_ws_root(json_cfg.is_ssd), json_cfg.test_name)
cfg_json_root = os.path.join(ppenv.get_cfg_json_ws(), json_cfg.test_name)

callsets = os.path.join(ppenv.get_callsets_root(), json_cfg.callsets)
vid_mapping = os.path.join(ppenv.get_mapping_vid_root(), json_cfg.vid_mapping)
ref_geno = os.path.join(ppenv.get_reference_root(), json_cfg.ref_geno)
#columns_list = get_query_pos_list(input_fn, region_pos, json_cfg.array_start_pos)

for isample in range(3):
    nsk = num_samples[isample]           # num samples in k
    exp_json_cfg_root = os.path.join(cfg_json_root, 'prep4interestpos%d_%d' % (num_parallel, nsk))
    if not os.path.isdir(exp_json_cfg_root):
        os.makedirs(exp_json_cfg_root)
    n_name = "%s_prep%dk" % (json_cfg.test_name, nsk)
    ub = num_samples[isample] - 1
    num_array = len(json_cfg.array_start_pos)
    tdb_ws = os.path.join(tdb_ws_root, "%d_%d-%s" % (num_array, num_parallel, n_name))
    col_part = CREATE_COL_PARTITION_F(json_cfg.array_start_pos, tdb_ws)
    nt_loader_cfg_t = LOADER_CFG_T()
    nt_loader_cfg = nt_loader_cfg_t._replace(column_partitions=col_part, \
    callset_mapping_file=callsets, vid_mapping_file=vid_mapping, \
    reference_genome=ref_geno, delete_and_create_tiledb_array=False, \
    lb_callset_row_idx=0, ub_callset_row_idx=ub, size_per_column_partition=(ub+1)*__COL_PARTITION_SIZE_UNIT)

    json_ldfn = os.path.join(exp_json_cfg_root, "%d-lc_%s_prep.json" % (num_array, json_cfg.test_name))
    with open(json_ldfn, 'w') as ofd:
        json.dump(nt_loader_cfg._asdict(), ofd, sort_keys=True, indent=4)
    print("INFO: generated LOADER", json_ldfn)

    nt_fsquery_t = QUERY_FULL_SCAN_CFG_T(json_cfg.num_parallel)
    my_nt_query_fs = nt_fsquery_t._replace(workspace=tdb_ws, \
        query_attributes=IllmnQueryAttr, \
        reference_genome=ref_geno, \
        vid_mapping_file=vid_mapping, \
        callset_mapping_file=callsets)
    for i in range(json_cfg.num_parallel):
        my_nt_query = my_nt_query_fs._replace(array=nt_loader_cfg.column_partitions[i]["array"])
        json_qtfn = os.path.join(exp_json_cfg_root, "fsqc_%s_%d.json" % (json_cfg.test_name, i))
        with open(json_qtfn, 'w') as ofd:
            json.dump(my_nt_query._asdict(), ofd, sort_keys=True, indent=4)
        print("INFO: generated QUERIER", json_qtfn)
