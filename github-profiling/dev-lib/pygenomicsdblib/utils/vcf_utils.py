#! /usr/bin/python3
# pylint: disable=too-many-locals

import random
import os.path
from subprocess import Popen, call
from collections import OrderedDict
import json
import logging
import uuid
from shutil import copyfile
import vcf

from . import constant_def

logger = logging.getLogger()
dir_limit = 4294967295

def generate_histogram(exec_name, vid_fn, callsets_fn, hist_fn):
    nt_hist = constant_def.HISTOGRAM_CFG_T()
    hist_cfg = nt_hist._replace(callset_mapping_file=callsets_fn, vid_mapping_file=vid_fn)
    hist_cfg_fn = os.path.join(os.path.dirname(callsets_fn), os.path.basename(callsets_fn).replace("callsets", "hist"))
    with open(hist_cfg_fn, 'w') as ofd:
        json.dump(hist_cfg._asdict(), ofd, sort_keys=True, indent=4)
    logger.debug("generate_histogram made json config  %s", hist_cfg_fn)
    with open(hist_fn, 'w') as ofd:
        call([exec_name, hist_cfg_fn], shell=False, stdout=ofd)
    logger.info("generate_histogram made histogram %s", hist_fn)

def generate_vcfs(dest_dir, num_gen_files, from_vcf_files):
    ''' TODO: No longer in use. from_vcf_files is a list of vcf file path '''
    origial_files = [f for f in from_vcf_files if os.path.isfile(f) and not f.startswith(".claimed.")]
    input_len = len(origial_files)
    dest_file_prefix = os.path.join(dest_dir, 'gen%d' % input_len)
    n_batch = 1
    file_map = []
    sidx = len([os.path.join(dest_dir, fn) for fn in os.listdir(dest_dir) if fn[-3:] == ".gz"])
    while num_gen_files > 0:
        num_files = min(num_gen_files, dir_limit)
        # dest_dir = os.path.join(dest_dir, "vcfs_%d" % n_batch)
        dest_bi = num_files*(n_batch-1) + sidx
        file_map.extend([__copy_file(dest_file_prefix, i + dest_bi, origial_files, random.randint(0, input_len-1)) for i in range(num_files)])
        num_gen_files -= num_files
        n_batch += 1
    logger.debug("file_map=%s", file_map[:5])
    return file_map

def __copy_file(dest_file_prefix, dest_idx, origial_files, src_idx):
    assert src_idx < len(origial_files), "index %d is >= len(origial_files) %d" % (src_idx, len(origial_files))
    src_fn = origial_files[src_idx]
    dest_fn = "%s_%d%s" % (dest_file_prefix, dest_idx, src_fn[-3:])
    copyfile(src_fn, dest_fn)
    # call(['rsync', '-a', src_fn, dest_fn])
    src_tbi_fn = "%s.tbi" % src_fn
    if os.path.isfile(src_tbi_fn):
        copyfile(src_tbi_fn, "%s.tbi" % dest_fn)
        #call(['rsync', '-a', src_tbi_fn, "%s.tbi" % dest_fn])
    return dest_fn, src_fn        

# callsets
def generate_callsets_json(vcf_inputfiles, callsets_fn):
    global_callset_idx = 0
    callsets_dict = OrderedDict()
    for vcf_file in vcf_inputfiles:
        read_mode = 'r' if vcf_file[-3:] == 'vcf' else 'rb'
        logger.debug("generate_callsets_json: processing %s", vcf_file)
        with open(vcf_file, read_mode) as fd:
            vcf_reader = vcf.Reader(fd)
            local_callset_idx = 0
            for callset_name in vcf_reader.samples:
                curr_callset_info = OrderedDict()
                if callset_name in callsets_dict:
                    ss = str(uuid.uuid4())
                    logger.warning('Duplicate callset name %s: appending _%s', callset_name, ss) 
                    callset_name += ('_' + ss)
                curr_callset_info["row_idx"] = global_callset_idx
                # curr_callset_info["idx_in_file"] = local_callset_idx
                curr_callset_info["filename"] = vcf_file
                callsets_dict[callset_name] = curr_callset_info
                local_callset_idx += 1
                global_callset_idx += 1
    with open(callsets_fn, 'w') as ofd:
        json.dump({'callsets' : callsets_dict, "file_division" : [vcf_inputfiles]}, ofd, indent=4, separators=(',', ': '))
    return global_callset_idx

# not tested for latest runs. may need adjustment
def append_callsets(origial_files, callsets_fn, gen_num):
    ''' update origial_files is a list of vcf file path '''
    assert os.path.isfile(callsets_fn)

    origial_files = [f for f in origial_files if os.path.isfile(f)]
    assert origial_files

    with open(callsets_fn, 'r') as ifd:
        callsets_dict_all = json.load(ifd)
    prev_callsets_dict = callsets_dict_all['callsets']
    start_idx = len(prev_callsets_dict)
    logger.info("gen_callsets_json_random found start_idx=%d, in callsets file %s", start_idx, callsets_fn)
   
    callsets_dict = OrderedDict(prev_callsets_dict)
    num_gen_i = start_idx
    curr_callset_info = OrderedDict()
    step = gen_num//10 if gen_num > 30 else 4
    for vcf_file in origial_files:
        read_mode = 'r' if vcf_file[-3:] == 'vcf' else 'rb'
        logger.debug("generate_callsets_json: processing %s", vcf_file)
        with open(vcf_file, read_mode) as fd:
            vcf_reader = vcf.Reader(fd)
            for local_idx, callset_name in enumerate(vcf_reader.samples):
                curr_callset_info = OrderedDict()
                if callset_name in callsets_dict:
                    ss = str(uuid.uuid4())
                    logger.warning('Duplicate callset name %s: appending _%s', callset_name, ss) 
                    callset_name += ('_' + ss)
                curr_callset_info["row_idx"] = num_gen_i
                curr_callset_info["idx_in_file"] = local_idx
                curr_callset_info["filename"] = vcf_file
                callsets_dict[callset_name] = curr_callset_info
                num_gen_i += 1
                if num_gen_i % step == 0:
                    last_2 = {k:callsets_dict[k] for i, k in enumerate(callsets_dict) if i >= ss and i < ss+2}
                    logger.debug("num_gen_i=%d, last_2=%s", num_gen_i, last_2)
        if num_gen_i - start_idx >= gen_num:
            break 
    with open(callsets_fn, 'w') as ofd:
        json.dump({'callsets' : callsets_dict}, ofd, indent=4, separators=(',', ': '))
    return num_gen_i
    
# launcher
def launch_vcf2tiledb(result_fn, lc_file_path):
    with open(result_fn, 'w') as fd:
        h_exec = Popen(['/usr/bin/time', '-v', 'vcf2tiledb', lc_file_path], shell=False, stdout=fd, stderr=fd)
        h_exec.communicate()
        return h_exec.wait()

# create callsets without copying vcf files to another folder. 
# NOT IN USE due to system caching
def create_callsets(dest_file_prefix, num_gen_files, from_vcf_files):
    ''' NOT IN USE from_vcf_files is a list of vcf file path '''
    origial_files = [f for f in from_vcf_files if os.path.isfile(f)]
    dest_dir = os.path.dirname(dest_file_prefix)
    #check_writable_dir(dest_dir)

    input_len = len(origial_files)
    dest_file_prefix = os.path.join(dest_dir, 'gen%d' % num_gen_files)
    n_batch = 1
    file_map = []
    while num_gen_files > 0:
        num_files = min(num_gen_files, dir_limit)
        dest_dir = os.path.join(dest_dir, "vcfs_%d" % n_batch)
        dest_bi = num_files*(n_batch-1)
        file_map.extend([__copy_file(dest_file_prefix, i + dest_bi, origial_files, random.randint(0, input_len-1)) for i in range(num_files)])
        num_gen_files -= num_files
        n_batch += 1
    logger.debug("file_map=%s", file_map[:5])
    return file_map
