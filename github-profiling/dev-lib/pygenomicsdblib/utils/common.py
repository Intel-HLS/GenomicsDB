#! /usr/bin/python3
#pylint: disable=broad-except

import logging
import shutil
import sys
import os
import os.path
import time
import platform
import inspect
from functools import wraps
from subprocess import call
import stat
from .constant_def import DEFAULT_WORKSPACE, LOGGER_ROOT

def extra_func_info(f):
    ''' wrap for print more info'''
    # see https://stackoverflow.com/questions/43506378/how-to-get-source-code-of-function-that-is-wrapped-by-a-decorator
    # for print print(inspect.getsource(f.__wrapped__)) - it print the source code
    @wraps(f)
    def wrapper(*args, **kwds):
        func_name = f.__name__
        launcher = os.path.basename(sys.modules['__main__'].__file__)
        caller_name = inspect.getmodule(inspect.stack()[1][0]).__name__
        caller_func = inspect.stack()[1][3]
        print("++%s() launcher=%s, caller=%s:%s()" % (func_name, launcher, caller_name, caller_func))
        return f(*args, **kwds)
    return wrapper
   
def remake_path(target_path):
    if os.path.exists(target_path):
        shutil.rmtree(target_path)
    os.makedirs(target_path)

def get_stats_path(root_path):
    return os.path.join(root_path, 'stats_row'), os.path.join(root_path, 'stats')

def str2num(x):
    ''' bery forgiven '''
    try:
        return int(x)
    except Exception:
        try:
            return float(x)
        except Exception:
            return None

def get_default_workdir():
    working_dir = os.environ.get('WS_HOME', DEFAULT_WORKSPACE)
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
    return working_dir
    
def gen_timestamp():
    return int(time.time())

def check_writable_dir(dir_path, create_if_not=True):
    if not os.path.isdir(dir_path) and create_if_not:
        os.makedirs(dir_path)
        call(['chmod', '-R', '+w', dir_path])
        
def is_windows():
    return platform.system() == 'Windows'

def check_x_mode(file_path):
    if not is_windows():
        assert os.path.isfile(file_path), ' file %s not found' % file_path
        st = os.stat(file_path)
        os.chmod(file_path, st.st_mode | stat.S_IEXEC)

def get_my_logger(log_name, log_level=logging.INFO, log_to=2):
    '''log_to = 0, console, 1 - lof file, 2 - both '''
    format_str = '%(asctime)s %(name)s::%(levelname)s: %(message)s'
    handlers = []
    if log_to != 1:
        handlers.append(logging.StreamHandler(sys.stdout))
        print('INFO: log to console')
    if log_to != 0:
        fn = os.path.join(LOGGER_ROOT, "%s.%s.log" % (log_name, platform.node().split('.')[0]))
        if os.path.isfile(fn):
            os.remove(fn)
        handlers.append(logging.FileHandler(filename=fn))
        print('INFO: log file name is', fn)
    logging.basicConfig(format=format_str, datefmt='%m/%d %H:%M:%S', level=log_level, handlers=handlers)
    return logging.getLogger(log_name)        

# def get_console_logger(log_name, log_level=logging.INFO):
#     format_str = '%(asctime)s %(name)s::%(levelname)s: %(message)s'
#     logging.basicConfig(format=format_str, datefmt='%y%m%d %H:%M:%S', level=log_level)
#     return logging.getLogger(log_name)

## convert iostat log into csv
def __write2file(cvs_prefix, pid, header, lines):
    ofile = "%s_%s.cvs" % (cvs_prefix, pid)
    with open(ofile, 'w') as ofd:
        ofd.write("%s\n" % ','.join(header))
        for fields in lines:
            ofd.write("%s\n" % ','.join(fields))
    return ofile

def make_csv_pidstat(from_log, to_csv_prefix):
    ''' head of all fields:
    Time, UID, PID, %usr, %system,  %guest, %CPU, CPU, minflt/s, majflt/s. VSZ, RSS, %MEM, kB_rd/s, kB_wr/s kB_ccwr/s, cswch/s, nvcswch/s, Command '''
    with open(from_log, 'r') as fd:
        lines = fd.readlines()
    f_idx = [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17]
    header = [f for i, f in enumerate(lines[2][1:].split()) if i in f_idx]
    dataline = [l.split() for l in lines[3:] if l[0] != '#' and len(l) > 20]
    pid_d = {}
    for sl in dataline:
        k = (sl[1], sl[2])
        v = [sl[i] for i in f_idx]
        if k in pid_d:
            pid_d[k].append(v)
        else:
            pid_d[k] = [v]
    pid_csv_files = [__write2file(to_csv_prefix, key[1], header, val) for  key, val in pid_d.items()]
    return pid_csv_files

def make_csv_iostat(from_log, to_csv_prefix):
    ''' head of all fields:
    Device:, rrqm/s, wrqm/s, r/s, w/s, rkB/s, wkB/s, avgrq-sz, avgqu-sz, await, r_await, w_await, svctm, %util'''
    csvfile = "%s.cvs" % (to_csv_prefix)
    f_idx = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13]    # skipped svctm
    has_header = False
    def write_out(ofd, line):
        data = [f for i, f in enumerate(line.split()) if i in f_idx]
        ofd.write("%s\n" % ','.join(data))

    with open(csvfile, 'w') as ofd:
        with open(from_log, 'r') as fd:
            for line in fd:
                if line.startswith('Device:'):
                    if not has_header: 
                        write_out(ofd, line)
                        has_header = True
                elif has_header:
                    write_out(ofd, line)
    return csvfile
