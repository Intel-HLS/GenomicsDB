#! /usr/bin/python3
import os
import os.path
from subprocess import check_output
from time import sleep, strftime, localtime
import platform
import argparse
import smtplib
import traceback
from email.mime.text import MIMEText
from psutil import process_iter, NoSuchProcess

#TDB_WS = "/mnt/app_ssd/scratch/illmn/16-x1_10000_illmn_varonly"

SLEEP_TIME = 30                     # 30 sec
MAX_NO_TARGET = 4                   # 2 min 
PROC_NAMES = ['gt_mpi_gather', 'vcf2tiledb']
LOG_ROOT = '/path/to/proc_logs'

class ProcessInfoCollector():
    PROC_ATTR = ['pid', 'status', 'create_time', 'name', 'cpu_percent', 'memory_percent', 'num_fds', 'num_threads', 'cpu_num', 'io_counters', \
    'threads', 'terminal', 'gids', 'memory_info', 'exe', 'uids', 'username', 'cwd', 'cmdline', \
    'connections', 'ppid', 'num_ctx_switches', 'open_files']
    receiver = 'mingx.rutar@intel.com'

    def __init__(self, fname, interval, tdb_ws):
        self.attr_div = self.PROC_ATTR.index('io_counters')
        self.sleep_time = interval
        self.node_id = platform.node().split('.')[0].split('-')[-1]
        self.sender = 'profiler-%s@sparkdmz.com' % self.node_id
        self.tdb_ws = tdb_ws
        ts = strftime("%y%m%d%H%M", localtime())
        self.fn = fname if fname else os.path.join(LOG_ROOT, 'pinfo_%s-%s.log' % (self.node_id, ts))
        with open(self.fn, 'w') as ofd:
            ofd.write('Start at %s, interval=%d, log to %s' % (ts, self.sleep_time, self.fn))
        self.cNotFound = 0

    def log_proc_info(self, procs=PROC_NAMES, bprint=False, bprntmore=False):
        bFound = False
        hasLoader = False
        with open(self.fn, 'a') as ofd:
            timestr = strftime("%y%m%d %H:%M:%S", localtime())
            ofd.write('PPP ---- proc_info @ %s: %s \n' % (timestr, ','.join(procs)))
            print('---- proc_info @ %s: %s' % (timestr, ','.join(procs)))
            for proc in process_iter():
                try:
                    pinfo = proc.as_dict(attrs=self.PROC_ATTR)
                    if pinfo['name'] in procs:    # ['gt_mpi_gather', 'vcf2tiledb', 'mpirun', 'time', 'python3']:
                        bFound = True
                        open_files = pinfo.pop('open_files')
                        info_l1 = ["%s: %s" % (k, pinfo[k]) for k in self.PROC_ATTR[:self.attr_div]]
                        ofd.write("%s\n" % str(info_l1))
                        info_l2 = ["%s: %s" % (k, pinfo[k]) for k in self.PROC_ATTR[self.attr_div:len(pinfo)]]
                        ofd.write("%s\n" % str(info_l2))
                        # logfd.write('myproc_info: ---- %s(%s) ---- \n' % (pinfo['name'], pinfo['pid']))
                        if bprint:
                            print(info_l1)
                            if bprntmore:
                                print(info_l2)
                            # print('ofiles=%s' % ','.join(ofiles[:2]))
                        if pinfo['name'] == 'vcf2tiledb':
                            hasLoader = True
                except NoSuchProcess:
                    pass
                except Exception as ex:
                    ofd.write(str(ex))
                    ofd.write(traceback.format_exc())
                ofd.flush()
        self.cNotFound = 0 if bFound else 1 + self.cNotFound
        return hasLoader

    def log_command(self, cmdline, bprint=False):
        ''' run the command and log result
        '''
        with open(self.fn, 'a') as ofd:
            timestr = strftime("%y%m%d %H:%M:%S", localtime())
            ofd.write("CCC ++++ %s @ %s\n" % (cmdline, timestr))
            retval = check_output(cmdline.split(), shell=False, timeout=self.sleep_time).decode("utf-8")
            ofd.write(retval)
            if bprint:
                print("CCC ++++ %s @ %s" % (cmdline, timestr))
                print(retval)

    def run(self, my_procs):
        ''' loop until 5 times could not find any target
        '''
        print_time = round(5 * 60 / self.sleep_time)     # about every 5 min
        detail_time = 6 * print_time                     # about 30 min
        num_run = 0
        lastHasLoader = False
        self.cNotFound = 0
        while True:
            bprnt = num_run % (print_time) == 0            
            bdetail = num_run % (detail_time) == 3
            hasloader = self.log_proc_info(my_procs, bprnt, bdetail)
            if self.tdb_ws and (lastHasLoader or hasloader):
                mycmd = 'du -h %s' % self.tdb_ws if bdetail  else 'du -sh %s' % self.tdb_ws
                pinfo_log.log_command(mycmd, bprnt)
            lastHasLoader = hasloader
            num_run += 1
            if self.cNotFound > MAX_NO_TARGET:
                print("cNotFound = %d, get out" % self.cNotFound)
                break
            sleep(self.sleep_time)
        return num_run

    def send_notice(self, note):
        msg = MIMEText(note)
        msg['Subject'] = 'Done profiling @ %s' % self.node_id
        msg['From'] = self.sender
        msg['To'] = self.receiver
        smtps = smtplib.SMTP('localhost')
        smtps.sendmail(self.sender, [self.receiver], msg.as_string())
        smtps.quit()

if __name__ == '__main__':
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-f', type=str, help='file name')
    my_parser.add_argument('-i', type=int, default=SLEEP_TIME, help='interval in sec')
    my_parser.add_argument('-t', type=str, help='tiledb working directory')
    myargs = my_parser.parse_args()

    start_time = strftime("%y%m%d%H%M", localtime())
    target_proc = ['mpirun', 'vcf2tiledb', 'gt_mpirun']
    print('Start at {}, log to {}, check {} every {} sec'.format(start_time, myargs.f, ','.join(target_proc), myargs.i))
    pinfo_log = ProcessInfoCollector(myargs.f, myargs.i, myargs.t)
    n_loop = pinfo_log.run(target_proc)
    end_time = strftime("%y%m%d%H%M", localtime())
    msg_text = "At {} finished checking {} that started at {}, log to {}, ".format(end_time, ','.join(target_proc), start_time, pinfo_log.fn)
    print(msg_text)
    pinfo_log.send_notice(msg_text)
