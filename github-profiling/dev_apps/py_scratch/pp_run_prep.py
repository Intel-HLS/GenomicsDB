#! /usr/bin/python3
import os
import os.path
import platform
from shutil import copyfile
from subprocess import Popen, call, DEVNULL, CalledProcessError, check_output
import time
n_parallel = 16

#node_cfg = {'27' : 5, '28' : 10, '29' : 15}
node_cfg = {'26' : 1000, '27': 3000, '28' : 7500, '29' : 15000}
node = platform.node().split('.')[0][-2:]

nn = node_cfg[node]
tdb_ws_empty_file = "/path/to/__tiledb_workspace.tdb"
tdb_ws = "/path/to/16_16-illmn_prep%d" % (nn)
json_root = "/path/to/prep4interestpos%d_%d"  % (n_parallel, nn)
interesting_pos_toproot = "/path/to//interesting_pos"

loader_exec = 'vcf2tiledb'
def loadtdb(mytdbws):
    if not os.path.isdir(mytdbws):
        os.makedirs(mytdbws)
    dest_file = os.path.join(mytdbws, os.path.basename(tdb_ws_empty_file))
    copyfile(tdb_ws_empty_file, dest_file)

    l_cfg = os.path.join(json_root, "16-lc_illmn_prep.json")
    assert os.path.isfile(l_cfg)
    cmd = "mpirun -np %d %s %s" % (n_parallel, loader_exec, l_cfg)
    print('launching cmd=', cmd)
    rc = call(cmd.split(), stdout=DEVNULL, stderr=DEVNULL, shell=False)
    # pid = Popen(cmd.split(), stdout=DEVNULL, stderr=DEVNULL, shell=False)
    print('Done loadtdb rc=', rc)

def wait_4_completion(process_path, check_interval_min):
    proc_name = os.path.basename(process_path)
    SleepTime = check_interval_min * 60
    num_wake_up = 0
    while True:
        try:
            procs = check_output(['pidof', proc_name]).split()
            time.sleep(SleepTime)
            num_wake_up += 1
            print("INFO %s running for %s minutes, pid_list=%s" % (proc_name, check_interval_min * num_wake_up, procs))
        except CalledProcessError:
            print("!! proc_name %s  no longer run" % proc_name)
            break
    return

query_exec = '/home/mingrutar/docker-geno/Release/bin/gt_mpi_gather'
def querydb(mytdbws, output_root):
    if not os.path.isdir(mytdbws):
        raise FileExistsError('tiledb workspace %s not exist', mytdbws)
    if not os.path.isdir(output_root):
        os.makedirs(output_root)
    pid_list = {}
    for i in range(n_parallel):
        ofd = None
        try:
            q_cfg = os.path.join(json_root, "fsqc_illmn_%d.json" % i)
            output = os.path.join(output_root, 'interesting_pos_%d' % i)
            qcmd = [query_exec, '-j', q_cfg, '--produce-interesting-positions']
            print("qcmd=%s, output=%s" % (qcmd, output))
            ofd = open(output, 'a')
            pid = Popen(qcmd, stdout=ofd, stderr=ofd, shell=False, close_fds=True)
            pid_list[i] = pid
            print("completed %s: i=%d, output=%s" % (pid, i, output))
        except Exception as ex:
            pid_list[i] = None
            print("WARN caught exception:: %s." %  ex)   #sys.exc_info()
        finally:
            if ofd:
                ofd.close()

if __name__ == '__main__':
    print('tdb_ws=', tdb_ws)
    if not os.path.isdir(tdb_ws):
        loadtdb(tdb_ws)
    # wait_4_completion(loader_exec, 15)
    interestpos_root = os.path.join(interesting_pos_toproot, "sample%d_%s" % (nn, node))
    querydb(tdb_ws, interestpos_root)
    print('Done')
