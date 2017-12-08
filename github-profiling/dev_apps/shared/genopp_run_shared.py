#! /usr/bin/python3
# pylint: disable=too-many-locals, too-many-arguments, attribute-defined-outside-init, global-statement
# module name: genopp_run_shared
# SharedTestRun
import sys
import os
import os.path
import platform
from subprocess import call, Popen
import uuid
from time import localtime, strftime, sleep, time
import traceback
import logging
import json
from shutil import copyfile, rmtree
from collections import OrderedDict
from psutil import disk_partitions, process_iter, NoSuchProcess
from pygenomicsdblib.remote.stat_capturer import StatCapturer
from pygenomicsdblib.utils.constant_def import MPIRUN_PATH, CLEAN_CACHE_SCRIPT_PATH, L_MARKER, Q_MARKER
logger = logging.getLogger()

assert os.path.isfile(MPIRUN_PATH), "Cannot find %s" % MPIRUN_PATH
tdb_ws_empty_file = os.path.join(os.path.dirname(__file__), "__tiledb_workspace.tdb")
assert os.path.isfile(tdb_ws_empty_file)
VCF_exec_root = os.path.join(os.path.expanduser('~'), 'docker-geno', 'Release', 'bin')

#  pp_run_prep.py uses '--produce-interesting-positions' & t_illmn_prep.py gen *[lq]c*.json
query_exec, loader_exec, result_log, output_fullpath = None, None, None, None

SAMPLING_INTERVAL = 0                    # 1 sec, <=0 do not sample
DEFAULT_Q_SEG_SIZE = 10000               # 10K is default segment size
DevNull = open(os.devnull, 'wb', 0)
# sampling
__sampler = None                         # StatCapturer 'iostat'

# No longer needed, used for verifying output of --print-calls. Only verify position 1 and 10
__record_output = lambda qc: qc.split('_')[-2] in ['1', '10']
__output_fn = lambda qc: os.path.join(output_fullpath, os.path.basename(qc).replace('json', 'output'))

nodename = platform.node().split('.')[0]

def can_exec_geno():
    ''' hacking '''
    hostname = platform.node().split('.')[0]
    return hostname[:-1] == 'compute-2-2' or hostname[:-1] == 'dhls-skx-'

def __prepare_tdb_ws(tdb_ws):
    if not os.path.isdir(tdb_ws):
        os.makedirs(tdb_ws)
        dest_file = os.path.join(tdb_ws, os.path.basename(tdb_ws_empty_file))
        copyfile(tdb_ws_empty_file, dest_file)

__has_sudo_priv = True          # hacking!
def __run_clear_cache(logrd):
    global __has_sudo_priv
    def no_priv_run():
        logrd.write("!!WARN: no sudo priv, run %s without sudo\n" % CLEAN_CACHE_SCRIPT_PATH)
        print("!!WARN: no sudo priv, run %s without sudo"  % CLEAN_CACHE_SCRIPT_PATH)
        call([CLEAN_CACHE_SCRIPT_PATH], stdout=logrd, stderr=logrd, shell=False)
    if __has_sudo_priv:
        try:
            call(['/usr/bin/sudo', CLEAN_CACHE_SCRIPT_PATH], stdout=logrd, stderr=logrd, shell=False)
        except KeyboardInterrupt:
            __has_sudo_priv = False
            no_priv_run()
    else:
        no_priv_run()

MAX_PROC_RUN_TIME = 4 * 60 * 60     # no proc runs over 4 hours. 

class ProcessInfoCollector():
    PROC_ATTR = ['pid', 'status', 'cwd', 'create_time', 'name', 'cmdline', 'cpu_percent', 'memory_percent', 'num_fds', \
    'num_threads', 'open_files', 'terminal', 'gids', 'cpu_num', 'memory_info', 'exe', 'threads', 'uids', 'username', \
    'connections', 'ppid', 'num_ctx_switches', 'io_counters']
    PROC_NAMES = ['gt_mpi_gather', 'vcf2tiledb', 'mpirun']  

    def __init__(self, root_dir):
        timestr = strftime("%y%m%d%H%M", localtime())
        self.fn = os.path.join(root_dir, 'proc_info_%s.log' % timestr)
        with open(self.fn, 'w') as ofd:
            ofd.write('Start at %s' % timestr)

    def log_proc_info(self):
        with open(self.fn, 'a') as ofd:
            timestr = strftime("%y%m%d %H:%M:%S", localtime())
            ofd.write('---- proc_info @ %s: %s \n' % (timestr, ','.join(self.PROC_NAMES)))
            for proc in process_iter():
                try:
                    pinfo = proc.as_dict(attrs=self.PROC_ATTR)
                    if pinfo['name'] in self.PROC_NAMES:    # ['gt_mpi_gather', 'vcf2tiledb', 'mpirun', 'time', 'python3']:
                        ofd.write("%s\n" % str(pinfo))
                        # logfd.write('myproc_info: ---- %s(%s) ---- \n' % (pinfo['name'], pinfo['pid']))
                        # logfd.write("%s\n" % str(pinfo))
                except NoSuchProcess:
                    pass
                except Exception as ex:
                    ofd.write(str(ex))
                    ofd.write(traceback.format_exc())
                    ofd.flush()

def __myproc_info(logfd):
    proc_attr = ['name', 'pid', 'status', 'cwd', 'create_time', 'num_fds', 'cpu_num']
    for proc in process_iter():
        try:
            proc_info = proc.as_dict(attrs=proc_attr)
            if proc_info['name'] in ['mpirun', 'time', 'python3']:
                logfd.write("%s\n" % ["%s=%s" % (x, proc_info[x]) for x in proc_attr])
                # logfd.write('myproc_info: ---- %s(%s) ---- \n' % (pinfo['name'], pinfo['pid']))
                # logfd.write("%s\n" % str(pinfo))
        except NoSuchProcess:
            pass
        except Exception as ex:
            logfd.write(str(ex))
            logfd.write(traceback.format_exc())
            logfd.flush()

def __check_procs(procs, logfd):
    finished = []
    for np, ppinfo in procs.items():
        rc = ppinfo[0].poll()
        if rc != None:
            with open(ppinfo[2], 'r') as tmpfd:
                logfd.write(tmpfd.read())
            timestr = strftime("%y%m%d %H:%M:%S", localtime())
            logfd.write("\nrun_the_cmd @ %s: Done pid=%s np=%d, rc=%d, tmpfn=%s\n" % (timestr, ppinfo[0].pid, np, rc, ppinfo[2]))
            # print("run_the_cmd@%s: Done pid=%s np=%d of %d, rc=%d\n" % (timestr, proc.pid, np, num_parallel, rc))
            os.remove(ppinfo[2])
            finished.append(np)
    _ = [procs.pop(k) for k in finished]

THREAD_AFFINITY = True

def execun_with_mpi(mycmd, fn):
    tmpfd = open(fn, 'wb')
    pexec = Popen(mycmd.split(), shell=False, stdout=DevNull, stderr=tmpfd, close_fds=True)
    print('execun_with_mpi, cmd=', mycmd)
    if not pexec:
        raise RuntimeError("launch test failed")
    return pexec

def __run_the_cmd(cmdline, poll_interval, lnum_parallel, num_parallel, logfd):
    def run_with_mpirun():
        tmpfn = '.%s_%d.%s' % (nodename[-2:], num_parallel, str(uuid.uuid4()))
        return {num_parallel : (execun_with_mpi(cmdline, tmpfn), cmdline, tmpfn)}
    def run_with_rank():
        ret = OrderedDict()
        for np in range(num_parallel):
            cmd = "%s -r %d" % (cmdline, np % lnum_parallel)
            tmpfn = '.%s_%d-%d.%s' % (nodename[-2:], num_parallel, np, str(uuid.uuid4()))
            ret[np] = execun_with_mpi(cmd, tmpfn)
        return ret
    procs = run_with_mpirun() if THREAD_AFFINITY else run_with_rank()
    pids = [(nn, pp[0].pid, pp[1]) for nn, pp in procs.items()]
    timestr = strftime("%y%m%d %H:%M:%S", localtime())
    # print("run_the_cmd @ %s: launched %d proc, pids=%s\n" % (timestr, num_parallel, pids))
    logfd.write("run_the_cmd @ %s: Launched %d proc, pids=%s\n" % (timestr, num_parallel, pids))
    # __myproc_info(logfd)
    logfd.flush()
    start_time = time()
    # ts1 = start_time
    while procs:
        sleep(poll_interval)            # 1 sec
        __check_procs(procs, logfd)
        if not procs:
            logfd.write("\nrun_the_cmd @ %s: Completed all %d processes\n" % (timestr, num_parallel))
            break
        delta_ts2 = time() - start_time
        if delta_ts2 >= MAX_PROC_RUN_TIME:
            logfd.write("\nrun_the_cmd @ %s: Error command '%s' with %d parallel has been running for %d s, exceeded allowed %d s. Incompleted proc are %s\n" % (timestr, cmdline, num_parallel, delta_ts2, MAX_PROC_RUN_TIME, [(np, ppinfo[0].pid) for np, ppinfo in procs.items()]))
            raise TimeoutError("command '%s' with %d parallel has been running for %d s, exceeded allowed %d s" % (cmdline, num_parallel, delta_ts2, MAX_PROC_RUN_TIME))
    logfd.flush()

def __run_query(query_cf, querier, lnum_parallel, qnum_parallel, tdb_ws):
    ''' tdb_ws is used for determining iostat sampling device'''
    def get_sampling_device(mypath):
        mnt_point = os.sep.join(os.path.abspath(mypath).split(os.sep)[:3])
        for p in disk_partitions():
            if p.mountpoint == mnt_point:
                return p.device
        return None
    def run_sampling(cmds):
        target_device = get_sampling_device(tdb_ws)
        if target_device:
            __sampler.start(base_fn.split('.')[0], target_device)
            __run_the_cmd(cmds, poll_interval, lnum_parallel, qnum_parallel, rfd)
            __sampler.stop()

    n_repeat = querier['repeat'] if 'repeat' in querier else 1
    seg_size_list = querier['seg_size'] if 'seg_size' in querier else [DEFAULT_Q_SEG_SIZE]
    base_fn = os.path.basename(query_cf)
    poll_interval = 1 if int(base_fn.split('_')[2]) <= 250 else 4
    print("**** RUN __run_query: cfg=%s, n_repeat=%d, produce_opt=%s, seg_size=%s" % (base_fn, n_repeat, querier['produce_opt'], seg_size_list))
    with open(result_log, 'a') as rfd:
        cmd_line_no_time = "%s -np %d %s -j %s %s" % (MPIRUN_PATH, qnum_parallel, query_exec, query_cf, querier['produce_opt']) if THREAD_AFFINITY else "%s -j %s %s" % (query_exec, query_cf, querier['produce_opt'])
        timestr = strftime("%y%m%d %H:%M:%S", localtime())
        rfd.write("\n+%s %s RUN __run_query: n_repeat=%d, produce_opt=%s, seg_size=%s\n" % (Q_MARKER, timestr, n_repeat, querier['produce_opt'], seg_size_list))
        rfd.flush()
        if __record_output(query_cf):        # record the output for verification. not in use
            with open(__output_fn(query_cf), 'w') as ofd:
                cmd_line_no_time_rec = "%s -s %d" % (cmd_line_no_time, DEFAULT_Q_SEG_SIZE)
                print("INFO: record outputs for %s" %  cmd_line_no_time_rec)
                rfd.write("INFO: record outputs for %s\n" %  cmd_line_no_time_rec)
                call(cmd_line_no_time_rec.split(), stdout=ofd, stderr=rfd, shell=False)
        for ss in seg_size_list:
            cmd_line = "/usr/bin/time -v %s -s %s" % (cmd_line_no_time, ss)   #segment size
            for i in range(n_repeat):
                rfd.write("--Start %d @ %s: %s\n" % (i, strftime("%y%m%d %H:%M:%S", localtime()), cmd_line))
                rfd.flush()
                __run_clear_cache(rfd)
                # if i == 0 and '_1000_' in base_fn and ss == seg_size_list[-1]:
                if __sampler and i == 0 and ss == seg_size_list[-1]:
                    run_sampling(cmd_line)
                else:
                    __run_the_cmd(cmd_line, poll_interval, lnum_parallel, qnum_parallel, rfd)
                rfd.write("--Done %d @ %s: %s\n\n" % (i, strftime("%y%m%d %H:%M:%S", localtime()), cmd_line))
                rfd.flush()
                print("--Done %d: %s" % (i, cmd_line))
        rfd.write("-%s %s END __run_query: n_repeat=%d, produce_opt=%s, seg_size=%s\n" % (Q_MARKER, strftime("%y%m%d %H:%M:%S", localtime()), n_repeat, querier['produce_opt'], seg_size_list))
        rfd.flush()

def __run_pair(loader_cf, query_cf, tdb_ws, qnum_parallel, load_options, query_options, isadd=False):
    ''' recreate_tdbws = True - if tdb_ws already loaded, reload
    '''
    with open(result_log, 'a') as rfd:
        task_name = "lc=%s, qc=%s, ws=%s, recreate_tdbws=%s" % (loader_cf, query_cf, tdb_ws, load_options['force_reload_first'])
        if not load_options['share'] and not isadd and os.path.isdir(tdb_ws):
            rmtree(tdb_ws)
        def need_load(loader_json, gdb_meta):
            ''' for now only check TEST0 '''
            with open(gdb_meta) as fd1:
                meta_data = json.load(fd1)
            with open(loader_json) as fd2:
                loader_cfg = json.load(fd2)
            rfd.write("loader_cfg['ub_callset_row_idx']=%d, meta_data['max_valid_row_idx_in_array']=%d\n" % (int(loader_cfg['ub_callset_row_idx']), int(meta_data['max_valid_row_idx_in_array'])))
            rfd.flush()
            return int(loader_cfg['ub_callset_row_idx']) > int(meta_data['max_valid_row_idx_in_array'])
        meta_fpath = os.path.join(tdb_ws, 'TEST0', "genomicsdb_meta.json")
        to_load = need_load(loader_cf, meta_fpath) if os.path.isfile(meta_fpath) else True
        if to_load:    # not os.path.isdir(tdb_ws) or isadd:
            print("+++ RUN __run_pair.loader: %s, isadd=%s" % (task_name, isadd))
            rfd.write("\n+%s %s RUN __run_pair.loader: %s" % (L_MARKER, strftime("%y%m%d %H:%M:%S", localtime()), task_name))
            rfd.write("--cmd: /usr/bin/time', '-v', %s, %s\n" % (loader_exec, loader_cf))
            rfd.flush()
            __prepare_tdb_ws(tdb_ws)
            __run_clear_cache(rfd)
            call(['/usr/bin/time', '-v', MPIRUN_PATH, '-np', str(load_options['num_parallel']), loader_exec, loader_cf], stdout=rfd, stderr=rfd, shell=False)
            rfd.write("-%s %s END RUN __run_pair.loader: %s\n" % (L_MARKER, strftime("%y%m%d %H:%M:%S", localtime()), task_name))
            rfd.flush()
            print(" --du: du -h %s" % (tdb_ws))
            rfd.write(" DDD: du -h %s\n" % (tdb_ws))
            call(['du', '-h', tdb_ws], stdout=rfd, stderr=rfd, shell=False)
        else:
            print("+++ SKIP __run_pair.loading for %s" % (task_name))
            rfd.write("\n+++ SKIP __run_pair.loading for %s\n" % (task_name))
        rfd.flush()
    if query_cf:
        __run_query(query_cf, query_options, load_options['num_parallel'], qnum_parallel, tdb_ws)
    print(" DONE ... task ", task_name)

def __get_exec(file_name):
    exec_fn = os.path.join(VCF_exec_root, file_name)
    assert os.path.isfile(exec_fn) and os.access(exec_fn, os.X_OK)
    return exec_fn

def __print_exec_info(fd):
    fd.write(" vcf2tiledb version is ")
    fd.flush()
    call([loader_exec, '--version'], stdout=fd, stderr=fd, shell=False)
    fd.write(" gt_mpi_gather version is ")
    fd.flush()
    call([query_exec, '--version'], stdout=fd, stderr=fd, shell=False)
    fd.write(" ldd gt_mpi_gather: ")
    fd.flush()
    call(['ldd', query_exec], stdout=fd, stderr=fd, shell=False)

def __run_first_pair(full_cfg_fn, run_options, log_root):
    if not can_exec_geno():
        print("INFO: __run_first_pair cannot exec: full_cfg_fn=", full_cfg_fn, ', run_options=', run_options)
        return

    global query_exec, loader_exec, result_log, output_fullpath, __sampler
    query_exec = __get_exec('gt_mpi_gather')
    loader_exec = __get_exec('vcf2tiledb')

    json_fn = full_cfg_fn[0]
    ws_dir = os.path.dirname(json_fn)
    assert os.path.isdir(ws_dir), "wrong loader file %s" % json_fn

    timestr = strftime("%y%m%d%H", localtime())
    hostname = platform.node().split('.')[0]
    if not os.path.isdir(log_root):
        os.makedirs(log_root)
    result_log = os.path.join(log_root, "%s_%s.result" % (hostname, timestr))

    output_fullpath = os.path.join(log_root, 'query_output')
    if os.path.isdir(output_fullpath):
        rmtree(output_fullpath)
    os.makedirs(output_fullpath)

    sampling_root = os.path.join(log_root, 'sampling')
    if os.path.isdir(sampling_root):
        rmtree(sampling_root)
    os.makedirs(sampling_root)
    __sampler = StatCapturer(sampling_root, sampling_root, 'iostat', SAMPLING_INTERVAL) if SAMPLING_INTERVAL > 0 else None

    with open(result_log, 'w') as fd:
        tdb_ws = full_cfg_fn[3]
        fd.write("==== Start Test @ %s .. tiledb %s recreate_tdbws is %s, is tdb exist %s \n" % \
                (strftime("%y%m%d %H:%M:%S", localtime()), tdb_ws, run_options['loader']['force_reload_first'], os.path.isdir(tdb_ws)))
        __print_exec_info(fd)
        if run_options['loader']['force_reload_first'] and os.path.isdir(tdb_ws):
            rmtree(tdb_ws)
            fd.write("removed %s " % tdb_ws)
        fd.flush()

    monitor = os.path.join(os.environ['HOME'], "bin", "proc_monitor")
    if os.path.isfile(monitor):
        pi_fn = os.path.join(log_root, 'procinfo_%s-%s.log' % (hostname.split('-')[-1], timestr))
        with open(pi_fn, 'w') as ofd:
            ofd.write('Start at %s' % timestr)
        # pi_cmd = "/usr/bin/python3 %s -f %s" % (monitor, pi_fn)
        Popen(["/usr/bin/python3", monitor, "-f", pi_fn], shell=False)
    
    __run_config(full_cfg_fn[0], full_cfg_fn[1], full_cfg_fn[2], full_cfg_fn[3], full_cfg_fn[4], run_options['loader'], run_options['querier'])

def __run_config(base_fn, add_fn, query_fn, myws, qnum_parallel, load_options, query_options):
    with open(result_log, 'a') as rfd:
        if os.path.isfile(query_fn) and (base_fn or add_fn):
            if base_fn and os.path.isfile(base_fn):
                rfd.write("INFO: find loader config (%s) and query config (%s),...\n" % (base_fn, query_fn))
                __run_pair(base_fn, query_fn, myws, qnum_parallel, load_options, query_options)
            if add_fn and os.path.isfile(add_fn):
                __run_pair(add_fn, query_fn, myws, qnum_parallel, load_options, query_options, True)
                rfd.write("INFO: find add loader config (%s) and query config (%s),...\n" % (add_fn, query_fn))
        else:
            rfd.write("\nWARN: cannot find query config (%s), ignored \n" % (query_fn))
        rfd.flush()

def run_tests(config_files, run_options, log_root):
    ''' config_files - list of (loader_fn.json, loader_add_fn.json, query.json, tdb_ws)
        run_options = {'loader' : {'share' : True, 'force_reload_first' : False, 'num_parallel' : 16},
                'querier' : {'produce_opt' : '--print-calls', 'repeat' : 5}}
        'produce_opt' is required
        checkpoint tested_list is not implemented.
    '''
    assert config_files
    assert run_options and 'produce_opt' in run_options['querier'] and 'repeat' in run_options['querier']
    #and run_options['querier']['produce_opt'] in gt_mpi_gather_options

    first_cfg = config_files.pop(0)
    __run_first_pair(first_cfg, run_options, log_root)

    # run the rest
    for base_fn, add_fn, query_fn, myws, num_parallel in config_files:
        __run_config(base_fn, add_fn, query_fn, myws, num_parallel, run_options['loader'], run_options['querier'])
#        __run_config(base_fn, add_fn, query_fn, myws, num_parallel, not run_options['loader']['share'], run_options['querier'])

def dryrun_tests():
    global query_exec, loader_exec
    query_exec = __get_exec('gt_mpi_gather')
    loader_exec = __get_exec('vcf2tiledb')
    __print_exec_info(sys.stdout)

# def __run_the_cmd_orig(cmdline, num_parallel, logfd):
#     procs = OrderedDict()
#     for np in range(num_parallel):
#         cmd = "%s -r %d" % (cmdline, np)
#         tmpfn = '.%s_%d-%d.%s' % (nodename[-2:], num_parallel, np, str(uuid.uuid4()))
#         tmpfd = open(tmpfn, 'wb')
#         pexec = Popen(cmd.split(), shell=False, stdout=DevNull, stderr=tmpfd, close_fds=True)
#         if not pexec:
#             raise RuntimeError("launch test failed")
#         procs[np] = (pexec, cmd, tmpfn)
#     pids = [(nn, pp[0].pid, pp[1]) for nn, pp in procs.items()]
#     timestr = strftime("%y%m%d %H:%M:%S", localtime())
#     # print("run_the_cmd @ %s: launched %d proc, pids=%s\n" % (timestr, num_parallel, pids))
#     logfd.write("run_the_cmd @ %s: Launched %d proc, pids=%s\n" % (timestr, num_parallel, pids))
#     __myproc_info(PROC_NAMES[0] + PROC_NAMES[1], logfd)
#     logfd.flush()
#     for np, ppinfo in procs.items():
#         proc = ppinfo[0] 
#         rc = proc.wait()
#         if rc == 0:
#             with open(ppinfo[2], 'r') as tmpfd:
#                 logfd.write(tmpfd.read())
#             timestr = strftime("%y%m%d %H:%M:%S", localtime())
#             logfd.write("\nrun_the_cmd @ %s: Done pid=%s np=%d of %d, rc=%d, tmpfn=%s\n" % (timestr, proc.pid, np, num_parallel, rc, ppinfo[2]))
#             # print("run_the_cmd@%s: Done pid=%s np=%d of %d, rc=%d\n" % (timestr, proc.pid, np, num_parallel, rc))
#             os.remove(ppinfo[2])
#         else:
#             logfd.write("\nrun_the_cmd@%s: Error pid=%s np=%d of %d, rc=%d\n" % (timestr, proc.pid, np, num_parallel, rc))
#     logfd.flush()
