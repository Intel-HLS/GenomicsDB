#!/usr/bin/env python3
#pylint: disable=wrong-import-position

import os
from platform import system
from shutil import rmtree
import json
import argparse
import logging
from pygenomicsdblib.utils.common import get_my_logger
logger = get_my_logger(os.path.basename(__file__)[:-3], logging.DEBUG)

from pygenomicsdblib.utils.vcf_utils import generate_vcfs, generate_callsets_json, generate_histogram

# inputs
MINI_TEST = "test"
GENO1000 = "geno1000"
ILLMN = "illmn"
WIN_SET = "win-set"
WIN_COPY_SET = "win-copy-set"

_src_dst_dirs = {
    # value =  (src, dst_root)
    MINI_TEST: ("/path/to/vcfs", "/path/to/scratch"),
    GENO1000: ("/scratch/1000genome/vcfs/", "/scratch"),
    ILLMN: ("/data/scratch/kdatta1/illmn/", None),
    WIN_SET: ("/path/to/vcfs", None),
    WIN_COPY_SET: ("/path/to/vcfs", "/path/to/scratch_genvcfs")}

_dst_root_path = "path/to/dev_data"

vcf_hist_exec = "vcf_histogram"
illmn_vid_fn = "/path/to/illmn_vid.json"

def copy_vcfs(num_files, src_vcfs_dir, dst_vcfs_root, remove_old=False):
    dst_vcfs_dir = os.path.join(dst_vcfs_root, "gen%d" % num_files)
    if os.path.exists(dst_vcfs_dir):
        if not remove_old:
            num_dst_files = [os.path.join(dst_vcfs_dir, fn) for fn in os.listdir(dst_vcfs_dir) if fn[-3:] == ".gz"]
            num_dst_idx_files = [os.path.join(dst_vcfs_dir, fn) for fn in os.listdir(dst_vcfs_dir) if fn[-7:] == ".gz.tbi"]
            if num_dst_files == num_dst_idx_files and num_dst_idx_files == num_files:
                return dst_vcfs_dir
        rmtree(dst_vcfs_dir)
    os.makedirs(dst_vcfs_dir)

    src_vcf_files = [os.path.join(src_vcfs_dir, fn) for fn in os.listdir(src_vcfs_dir) if fn[-3:] == ".gz" and not fn.startswith(".claimed.")]
    file_map = generate_vcfs(dst_vcfs_dir, num_files, src_vcf_files)

    file_map_fn = "file_map_%d.json" % num_files
    with open(file_map_fn, "w") as fd:
        json.dump(file_map, fd)
    return dst_vcfs_dir

if __name__ == "__main__":
    ''' python3 t_gen_callsets_copy.py -s illmn -n 15000 --histogram hist_file_name -f
       -f is needed only if duplicate VCF files are needed, which is Illmn VCF generator is available '''
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-s', type=str, default='illmn', required=False, help='input source: geno1000 or illmn or test, default illmn')
    my_parser.add_argument('-n', type=int, default=1000, help='number vcf, if > #source, duplicate (full copy) vcf from source')
    my_parser.add_argument('--histogram', type=str, required=False, help='generate histograme')
    my_parser.add_argument('-f', type=bool, default=False, help='recopy if exist, only when -g >0, default is False')
    args = my_parser.parse_args()

    logger.debug('args are %s', args)
    src_dir, dst_root = _src_dst_dirs[WIN_SET if system() == 'Windows' else args.s]
    logger.debug('src_dir=%s, dst_root=%s', src_dir, dst_root)
    vcfs_files = [os.path.join(src_dir, fn) for fn in os.listdir(src_dir) if fn[-3:] == ".gz"]
    logger.debug('# vcfs_files=%d', len(vcfs_files))
    if len(vcfs_files) < args.n and dst_root:
        logger.info("only %d original vcfs available, will generate vcfs @ %s", len(vcfs_files), dst_root)
        gen_vcfs_dir = copy_vcfs(args.n, src_dir, dst_root, args.f)
        vcfs_files = [os.path.join(gen_vcfs_dir, fn) for fn in os.listdir(gen_vcfs_dir) if fn[-3:] == ".gz"]
    total_num_vcf = min(args.n, len(vcfs_files))
    callsets_dir = os.path.join(_dst_root_path, args.s)
    callsets_fname = os.path.join(callsets_dir, "callsets_%s_%d.json" % (args.s, total_num_vcf))
    # the vcfs_files likely is out of order. we pick number of files first, then sort
    input_list = sorted(vcfs_files[:total_num_vcf], key=lambda x: int(os.path.basename(x).split('.')[0]))
    logger.debug("generate callsets size %d, #input=%d, #file=%d", total_num_vcf, len(input_list), len(vcfs_files))
    num_callsets = generate_callsets_json(input_list, callsets_fname)
    if args.histogram and os.path.isfile(callsets_fname):
        histogram_fname = os.path.join(callsets_dir, args.histogram)
        logger.debug('creating histogram callsets_fname=%s, histogram_fname=%s...', callsets_fname, histogram_fname)
        generate_histogram(vcf_hist_exec, illmn_vid_fn, callsets_fname, histogram_fname)
    else:
        logger.info("did not generate histogram: %s", 'OK' if not args.histogram else 'NO callsets %s' % callsets_fname)
    logger.info("DONE created %s with %d callsets, num_sample=%d", callsets_fname, total_num_vcf, len(vcfs_files))
