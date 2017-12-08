#! /usr/bin/python3

import os
import os.path
import time
from subprocess import call
import utils.constant_def as const
from utils.common import check_x_mode, remake_path

# remove to disable mpirun checking
def check_libraries(libs=None):
    ''' list of full path of expected_paths '''
    expected_paths = libs + const.DEFAULT_LIB_PATHS if libs else const.DEFAULT_LIB_PATHS
    if const.ENV_LIB_PATH in os.environ:
        print("check_libraries-1: %s=%s, expected_paths=%s" % (const.ENV_LIB_PATH, os.environ[const.ENV_LIB_PATH], expected_paths))
        ld_libs = os.environ[const.ENV_LIB_PATH].split(':')
        libs = expected_paths + [l for l in ld_libs if l not in expected_paths]
        os.environ[const.ENV_LIB_PATH] = ":".join(libs)
    else:
        print("check_libraries-2: expected_paths=%s" % (expected_paths))
        os.environ[const.ENV_LIB_PATH] = ":".join(expected_paths)
    print("check_libraries-3: %s=%s" % (const.ENV_LIB_PATH, os.environ[const.ENV_LIB_PATH]))

exec_list = ["/usr/bin/time", const.MPIRUN_PATH, const.CLEAN_CACHE_SCRIPT_PATH]   
def check_executions(cmd):
    exec_list.append(cmd)
    for exec_cmd in exec_list:
        if not os.path.exists(exec_cmd):
            raise RuntimeError('required command %s not found' % exec_cmd)
    
    os.system('pkill ' + os.path.basename(cmd))

def check_tdb_ws(is_loading, tdb_ws):
    if is_loading:
        remake_path(tdb_ws)
        for fn in ['__tiledb_workspace.tdb', '__tiledb_group.tdb']:
            with open(os.path.join(tdb_ws, fn), 'w') as fd:
                pass
    elif not os.path.isdir(tdb_ws):         #tdb_ws must exist for query 
        raise RuntimeError("required tiledb ws %s not found" % tdb_ws)
    assert os.path.isfile(const.CLEAN_CACHE_SCRIPT_PATH), 'cannot find %s' % const.CLEAN_CACHE_SCRIPT_PATH
    retcode = call(['sudo', const.CLEAN_CACHE_SCRIPT_PATH])
    if retcode == 0:
        time.sleep(const.WAIT_TIME_FOR_LAUNCHING)
    else:
        raise RuntimeError('Failed to run %s' % const.CLEAN_CACHE_SCRIPT_PATH)
        