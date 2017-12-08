#! /usr/bin/python3
#! pylint: disable=too-many-locals
# module name : genopp_process_utils
'''
Input data 1000g_interesting_positions is a supplied csv file (by Karthik). 1000g_interesting_positions is processed from 1000Genomics VCF data. The csv file has 4 columns: position, num_intersect, num_ref_block_cell and num_spos_aligned
position - bp position
num_intersect - number of samples intesect at this position
num_ref_block_cell - numbet of intesected positions that have reference
num_spos_aligned - number of positions that are aligned at the position
The number of variant is num_variant = num_ref_block_cell - num_intersect. In the other words, the intersected positions not defined in references.
The input file 1000g_interesting_positions contains 48604666 positions. About 3% (1486310) of them has variants.

copied from jupyter notebook @ C:\\Users\\mrutarx\\vezba\\DensityVisualization.ipynb '''

import json
import os
import os.path
import time
import glob
import logging
from pygenomicsdblib.utils.constant_def import DENSITY_NAMES, QUERY_CFG_T, LOADER_CFG_T, __COL_PARTITION_SIZE_UNIT
import profiling_env as ppenv

logger = logging.getLogger()
# def get_variant_df(src_fn):   not needed
#     ''' if the dataframe is saved, load from checkpoint. otherwise load from file and save df to checkpoint
#         src_fn: input file name
#     '''
#     assert src_fn
#     ckpt_root = os.path.join(ppenv.get_cfg_json_ws(), "checkpoints")
#     if not os.path.isdir(ckpt_root):
#         os.makedirs(ckpt_root)
#     base_name = os.path.basename(src_fn).split('.')
#     base_name = base_name if len(base_name) == 1 else base_name[:-1]
#     checkpoint_fn = os.path.join(ckpt_root, ''.join(base_name) +'.ckpt')
#     df = None
#     if os.path.isfile(checkpoint_fn):
#         if os.path.getmtime(src_fn) < os.path.getmtime(checkpoint_fn):
#             df = pd.read_pickle(checkpoint_fn)
#             assert isinstance(df, pd.core.frame.DataFrame)
#             print("INFO: loaded from checkpoint_fn=%s, df_type=%s" % (checkpoint_fn, type(df)))
#         else:
#             print('INFO source file %s has changed' % src_fn)
#             os.remove(checkpoint_fn)
#     if df is None:
#         df = __load_from_file(src_fn, checkpoint_fn)
#     return df
def __clean_dir(dir_name, glob_fl):
    if os.path.isdir(dir_name):
        for fl in glob_fl:
            os.remove(fl)           #   shutil.rmtree(dir_name)
    else:
        os.makedirs(dir_name)

#def gen_cfg_json4density(columns_list, share_tdb, json_cfg):
def gen_cfg_json4density(exp_name, is_test, pos_selector, json_cfg, share_tdb=True):
    ''' generates loader and query configs for denseity tests
        columns_list - list of [([pos..], [#var]), ..,(.., ..)]
        share_tdb - use one tdb. all query use this tdb
        json_cfg = test_name, callsets, vid_mapping, ref_geno
      Returns run_list
        run_list: loader.json, add.json, query.json, tdb_ws, num_parallel
        json_fn= {num_parall}-{ls|qc}_pos_numOfPos_lowRange-highRange.json
        if loader is shared, {num_parall}-ls_pos_TEST_NAME.json
    '''
    print()
    ppenv.set_run_type(is_test)
    print("gen_cfg_json4density: is_test=%s ppenv.is_test()=%s" % (is_test, ppenv.is_test()))
    assert pos_selector and json_cfg
    num2gen = json_cfg.num_pos2gen
    test_name = json_cfg.test_name
    assert test_name and num2gen
    region_pos = [num2gen] * len(DENSITY_NAMES)
    print("INFO: gen_cfg_json4density: region_pos=%s, json_cfg=%s, share_tdb=%s, is_test=%s" % (region_pos, json_cfg, share_tdb, ppenv.is_test()))
    num_array = len(json_cfg.array_start_pos)

    def __creat_loader(n_name, lb, ub):
        ''' ub is 5000 not 4999
        '''
        tdb_ws = os.path.join(tdb_ws_root, "%d-%s_%s" % (num_array, n_name, json_cfg.db_type))
        col_part = [{"array": "TEST%d" % (i), "begin": begin, "workspace": tdb_ws} \
            for i, begin in enumerate(json_cfg.array_start_pos)]
        nt_loader_cfg = nt_loader._replace(column_partitions=col_part, \
        callset_mapping_file=callsets, vid_mapping_file=vid_mapping, \
        reference_genome=ref_geno, delete_and_create_tiledb_array=False, \
        lb_callset_row_idx=lb, ub_callset_row_idx=ub-1, size_per_column_partition=(ub-lb)*__COL_PARTITION_SIZE_UNIT)
        json_ldfn = os.path.join(exp_json_cfg_root, "%d-lc_%s_%d-%d.json" % (num_array, n_name, lb, ub))
        # json_ldfn = os.path.join(cfg_json_root, exp_name, "%d-lc_%s.json" % (num_array, n_name))
        print("__creat_loader json_ldfn=", json_ldfn)
        with open(json_ldfn, 'w') as ofd:
            json.dump(nt_loader_cfg._asdict(), ofd, sort_keys=True, indent=4)
        print("INFO: __creat_loader generated LOADER", json_ldfn)
        return json_ldfn, tdb_ws, nt_loader_cfg

    def __create_query(n_name, position_lists, tdb_ws, nt_loader_cfg):
        array_list = [col_part["array"] for col_part in nt_loader_cfg.column_partitions[:json_cfg.num_parallel]]
        if len(array_list) != len(position_lists):
            print("!!! __create_query n_name=%s, #array_list=%d, #position_lists=%d" % (n_name, len(array_list), len(position_lists)))
            print("   array_list=%s \n   position_lists=%s" % (array_list, position_lists))
            return None
        query_cfg_nt = nt_query._replace(workspace=tdb_ws, array=array_list, query_column_ranges=position_lists, reference_genome=nt_loader_cfg.reference_genome, vid_mapping_file=nt_loader_cfg.vid_mapping_file, callset_mapping_file=nt_loader_cfg.callset_mapping_file, query_attributes=json_cfg.query_attributes)

        json_qtfn = os.path.join(exp_json_cfg_root, "%d-qc_%s.json" % (json_cfg.num_parallel, n_name))
        # json_qtfn = os.path.join(cfg_json_root, exp_name, "%d-qc_%s.json" % (num_array, n_name))
        with open(json_qtfn, 'w') as ofd:
            json.dump(query_cfg_nt._asdict(), ofd, sort_keys=True, indent=4)
        print("INFO: __create_query generated QUERIER", json_qtfn)
        return json_qtfn

    nt_loader = LOADER_CFG_T()
    nt_query = QUERY_CFG_T(json_cfg.num_parallel)
    tdb_ws_root = os.path.join(ppenv.get_tdb_ws_root(json_cfg.device), test_name)
    cfg_json_root = os.path.join(ppenv.get_cfg_json_ws(), test_name)

    run_list = []
    callsets = os.path.join(ppenv.get_callsets_root(), json_cfg.callsets)
    vid_mapping = os.path.join(ppenv.get_mapping_vid_root(), json_cfg.vid_mapping)
    ref_geno = os.path.join(ppenv.get_reference_root(), json_cfg.ref_geno)

    exp_json_cfg_root = os.path.join(cfg_json_root, exp_name)
    __clean_dir(exp_json_cfg_root, glob.iglob(os.path.join(exp_json_cfg_root, "*.json")))
    
    print("!! exp_json_cfg_root=", exp_json_cfg_root, ", exists=", os.path.isdir(exp_json_cfg_root))
    if share_tdb:   # 1 loader, N query
        l_name = "x%d_%d_%s" % (json_cfg.num_frag[0]+1, json_cfg.num_frag[1], json_cfg.test_name)
        num_samples, lb = 0, 0
        all_sel_pos = []
        for _ in range(json_cfg.num_frag[0]+1):          # num_frag[0] is num of additional
            lb = num_samples
            num_samples += json_cfg.num_frag[1]
            json_ldfn, tdb_ws, nt_loader_cfg = __creat_loader(l_name, lb, num_samples)
            df, _ = pos_selector.get_intpos_df(num_samples)
            columns_list, sel_info = pos_selector.get_query_pos_list(df, region_pos, json_cfg.array_start_pos, json_cfg.num_parallel)
            for i, region_cols in enumerate(columns_list):
                for cols in region_cols:
                    n_name = "%s_%d_%d-%d-%d" % (DENSITY_NAMES[i], len(cols[0][0]), min(cols[1]), max(cols[2]), num_samples)
                    json_qtfn = __create_query(n_name, cols[0], tdb_ws, nt_loader_cfg)
                    run_list.append((json_ldfn if lb == 0 else None, None if lb == 0 else json_ldfn, json_qtfn, tdb_ws, json_cfg.num_parallel))
            all_sel_pos.append(sel_info)
    else:      # N loader N query
        raise NotImplementedError
        # _ = [[__gen_json_func(i, cols, __creat_loader, __create_query) for cols in region_cols] for i, region_cols in enumerate(columns_list)]
    return run_list, all_sel_pos

def launch_profiling(run_list_and_opt, host_list, exp_name):
    ''' launch test remote or locally according to env
        run_list_and_opt - pickled run_list file and run options
        host_list - host list
    '''
    assert os.path.isfile(run_list_and_opt)
    assert host_list
    # in cfg_json_root/
    for host in host_list:
        exec_cmd, isremote = ppenv.get_remote_run_cmd(host, run_list_and_opt)
        run_dir = os.path.dirname(run_list_and_opt)
        timestr = time.strftime("%y%m%d%H", time.localtime())
        subdir_name = exp_name if exp_name else "run%s" % timestr
        log_root = os.path.join(run_dir, subdir_name, host)
        exec_cmd = "%s %s" % (exec_cmd, log_root)
        if not isremote:
            print("**launch locally cmd=", exec_cmd)
            os.system(exec_cmd)
        else:
            __clean_dir(log_root, glob.iglob(os.path.join(log_root, "nohup.*")))
            log_out_fn = os.path.join(log_root, "nohup.out")
            log_err_fn = os.path.join(log_root, "nohup.err")
            nohup_exec_cmd = "nohup %s >%s 2>%s &" % (exec_cmd, log_out_fn, log_err_fn)
            print("**launch remotely cmd=", nohup_exec_cmd)
            os.system(nohup_exec_cmd)
