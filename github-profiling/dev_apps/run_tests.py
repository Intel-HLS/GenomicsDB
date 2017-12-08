#! /usr/bin/python3
#! pylint: disable=invalid-name, line-too-long, wrong-import-position
''' ro run:
1) USE /home/mingrutar/bin/t_run_remote.bash for having nohup
2) @ local,
    nohup /usr/bin/python3 ./run_local.py -r N 1>/path/to/nohup_log.out 2>/path/to/nohup_log.err &

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

node_name = platform.node().split('.')[0]
# or 2 if want both console and file output   
logger = get_my_logger("%s-%s"%(os.path.basename(__file__)[:-3], node_name), logging.DEBUG, 0) 

target = 'illmn'                                           #variant-density/
run_root = "/path/to/run_tests/"
log_root = "/path/to/run_logs/"

RUN_SETTING = 6      # key in run_list_cfg
# for e-n0: my_run_list = lambda hid: [json_cfgs[0], json_cfgs[1], json_cfgs[-1]] if hid in [ "7", "9"] else [json_cfgs[2], json_cfgs[-2]]
#lambda hid: ["run_list-i-10", "run_list-i-30"] if hid in ["7", "9"] else ["run_list-i-20", "run_list-i-40"]
run_list_cfg = {
    # illmn, repeat item 2. due to ssd limit, only run 1, 4 mpiruns
    's' : lambda hid: ["run_list-s-%d0" % (x) for x in range(1, 5)],
    # illmn, repeat item 2. due to ssd limit, only run 1, 4 mpiruns
    't' : lambda hid: ["run_list-t-%d0" % (x) for x in range(1, 5)]
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
    assert hostname in ["compute-2-2%d" % (x) for x in range(7, 9)]
    rl = run_list_cfg[args.r](hostname[-1])
    logger.info("Start run_local: -r=%s, -d=%s run_list: %s", args.r, args.d, rl)
    print("Start run_local: -r={}, -d={} run_list: {}".format(args.r, args.d, rl))
    for cfg in rl:
        run_fn = os.path.join(run_root, target, cfg)
        with open(run_fn, 'rb') as ifd:                 
            run_cmd_cfg = pickle.load(ifd)  
        exp_name_list = sorted(run_cmd_cfg.keys(), key=lambda x: int(x[2:4]))
        _ = [run_exp(run_cmd_cfg[en], en, os.path.join(log_root, target), args.d) for en in exp_name_list]
       
        for exp_name, cmds in run_cmd_cfg.items():
            run_exp(cmds, exp_name, os.path.join(log_root, target), args.d)
