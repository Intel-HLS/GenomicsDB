#! /usr/bin/python3
''' set host-level environ.
        'add_to_lib_path' => add libs to LD_LIBRARY_PATH
'''
import sys
import os
import os.path
import pickle
import platform
from pprint import pprint
from pygenomicsdblib.utils.constant_def import ENV_LIB_PATH, DEFAULT_LIB_PATHS

# root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# if not root_path in sys.path:
#     sys.path.append(root_path)
from shared.genopp_run_shared import run_tests, dryrun_tests

def run_exp(run_cmd, my_exp_name="no_name"):
    print("String profiling %s ..." % my_exp_name)
    pprint(run_cmd)
    if 'add_to_lib_path' in run_cmd['run_options']:
        os.environ[ENV_LIB_PATH] = ':'.join(run_cmd['run_options']['add_to_lib_path']  + DEFAULT_LIB_PATHS)
    print("INFO: run_pickle.py: os.environ[%s] = %s" % (ENV_LIB_PATH, os.environ[ENV_LIB_PATH]))
    # show exec info launched locally
    if sys.stdin.isatty():
        dryrun_tests()
        print("INFO: Done run_pickelr locally", my_exp_name)
    else:
        print("\n INFO run_pickelr @ {}: launched from remote".format(platform.node()))
        run_tests(run_cmd['configs'], run_cmd['run_options'], log_root)
        print("INFO: Done profiling", my_exp_name)

if __name__ == '__main__':
    # run_pickle arg_list_fn
    assert len(sys.argv) == 3
    assert os.path.isfile(sys.argv[1])
    log_root = sys.argv[2]
    #is_running_ssh = 'SSH_CLIENT' in os.environ or 'SSH_TTY' in os.environ
    with open(sys.argv[1], 'rb') as ifd:
        run_cmd_cfg = pickle.load(ifd)
    if isinstance(run_cmd_cfg, dict):
        for exp_name, cmds in run_cmd_cfg:
            run_exp(cmds, exp_name)
    else:
        run_exp(run_cmd_cfg)
