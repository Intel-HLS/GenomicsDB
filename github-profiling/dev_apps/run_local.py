#! /usr/bin/python3
#! pylint: disable=invalid-name, line-too-long, wrong-import-position
''' ro run:
1) USE /home/mingrutar/bin/t_run_remote.bash for having nohup
2) nohup python3 ./run_local.py -r 4 1>/path/to/nohup_log.out 2>/path/to/nohup_log.err &

set host-level environ.
    'add_to_lib_path' => add libs to LD_LIBRARY_PATH
'''
import os
import os.path
import pickle
import platform
from pprint import pprint
import logging
import argparse
from pygenomicsdblib.utils.constant_def import ENV_LIB_PATH, DEFAULT_LIB_PATHS
from pygenomicsdblib.utils.common import get_my_logger
from shared.genopp_run_shared import run_tests, dryrun_tests
import shared.profiling_env as ppenv

node_name = platform.node().split('.')[0]
# or 2 if want both console and file output   
logger = get_my_logger("%s-%s"%(os.path.basename(__file__)[:-3], node_name), logging.DEBUG, 0) 

DEFAULT_TARGET = 'illmn'                                           #variant-density/

RUN_SETTING = 6      # key in run_list_cfg
# for e-n0: my_run_list = lambda hid: [json_cfgs[0], json_cfgs[1], json_cfgs[-1]] if hid in [ "7", "9"] else [json_cfgs[2], json_cfgs[-2]]
#lambda hid: ["run_list-i-10", "run_list-i-30"] if hid in ["7", "9"] else ["run_list-i-20", "run_list-i-40"]
run_list_cfg = {
    # e-10, .. e-50, parallel 1, 8,16, 44 and 88
    1 : lambda hid: ["run_list-e-%d0" % (x) for x in [1, 2, 5]] if hid in ["7", "9"] else ["run_list-e-%d0" % (x) for x in [3, 4]],
    # illmn  i-10,.. i-50, for 1, 4, 8, 16 mpirun; for 5k, 10k and 15k (actually all 1K, but)
    2 : lambda hid: ["run_list-i-10", "run_list-i-30"] if hid in ["7", "9"] else ["run_list-i-20", "run_list-i-40"],
    # illmn, repeat item 2. due to ssd limit, only run 1, 4 mpiruns
    3 : lambda hid: ["run_list-i-%d0" % (x) for x in range(1, 3)],
    # illmn, repeat item 2. due to ssd limit, only run 1, 4 mpiruns
    4 : lambda hid: ["run_list-a-%d0" % (x) for x in range(1, 3)],
    # illmn, repeat item 2. due to ssd limit, only run 1, 4 mpiruns
    5 : lambda hid: ["run_list-x-%d0" % (x) for x in range(1, 5)],
    # illmn, repeat item 2. due to ssd limit, only run 1, 4 mpiruns
    6 : lambda hid: ["run_list-y-%d0" % (x) for x in range(1, 5)],
    # illmn, repeat item 2. due to ssd limit, only run 1, 4 mpiruns
    7 : lambda hid: ["run_list-c-%d0" % (x) for x in range(1, 5)],
    # illmn, parallel 32, 44, 88 SSD only
    8 : lambda hid: ["run_list-d-%d0" % (x) for x in range(1, 4)],
    # illmn, parallel 16 optane only, on skx-0
    9 : lambda hid: ["run_list-o-%d0" % (x) for x in range(1, 2)],
    # illmn, parallel 16 hdd on skx-2
    10 : lambda hid: ["run_list-p-%d0" % (x) for x in range(1, 2)],
    # small 5 or 20 pos, --print-call
    101 : lambda hid: ["run_list-s-%d0" % (x) for x in range(1, 2)],
    # small 5 or 20 pos, --print-AC
    102 : lambda hid: ["run_list-t-%d0" % (x) for x in range(1, 2)],
    # small 5 or 20 pos, --print-AC
    103 : lambda hid: ["run_list-u-%d0" % (x) for x in range(1, 2)]
}

def run_exp(run_cmd, my_exp_name, my_log_root, dryrun=False):
    logger.info("profiling %s ...dryrun=%s", my_exp_name, dryrun)
    pprint(run_cmd)
    lib_paths = run_cmd['run_options']['add_to_lib_path']  + DEFAULT_LIB_PATHS \
        if 'add_to_lib_path' in run_cmd['run_options'] else DEFAULT_LIB_PATHS
    os.environ[ENV_LIB_PATH] = ':'.join(lib_paths)
    # show exec info launched locally
    test_log_dir = os.path.join(my_log_root, my_exp_name, node_name)
    logger.info("run_local @ %s, log_root at %s", node_name, test_log_dir)
    print("run_local @ %s, log_root at %s" % (node_name, test_log_dir))
    if not dryrun:
        run_tests(run_cmd['configs'], run_cmd['run_options'], test_log_dir)
        logger.info("Done profiling %s", my_exp_name)
        print("Done profiling %s" % my_exp_name)
    else:
        dryrun_tests()
        print("Done dryrun %s" % my_exp_name)

if __name__ == '__main__':
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d', type=bool, nargs='?', const=True, default=False, help='dryrun')
    my_parser.add_argument('-r', type=int, nargs='?', const=RUN_SETTING, default=RUN_SETTING, help='run config number, default %d' % RUN_SETTING)
    args = my_parser.parse_args()

    hostname = platform.node().split('.')[0]
    rl = run_list_cfg[args.r](hostname[-1])
    target = 'illmn-test' if args.r > 100 else DEFAULT_TARGET 
    logger.info("Start run_local: -r=%s, -d=%s, target=%s, run_list: %s", args.r, args.d, target, rl)
    print("Start run_local: -r={}, -d={}, target={} run_list: {}".format(args.r, args.d, target, rl))
    json_config_dir = os.path.join(ppenv.get_cfg_json_ws(), DEFAULT_TARGET)
    log_dir = os.path.join(ppenv.get_runlog_root(), target)
    for cfg in rl:
        run_fn = os.path.join(json_config_dir, cfg)
        with open(run_fn, 'rb') as ifd:                 
            run_cmd_cfg = pickle.load(ifd)  
        exp_name_list = sorted(run_cmd_cfg.keys(), key=lambda x: int(x[2:4]))
        _ = [run_exp(run_cmd_cfg[en], en, log_dir, args.d) for en in exp_name_list]
