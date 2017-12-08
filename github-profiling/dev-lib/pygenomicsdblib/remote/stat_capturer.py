#! /usr/bin/python3
# pylint: disable=broad-except

import os
import os.path
import time
from subprocess import Popen
import logging
# from ..utils import common

logger = logging.getLogger()

class StatCapturer():
    #  pidstat -hdIruw  -C vcf2tiledb 1
    stat_args = {'pidstat' : lambda exec_name: (0, ['/usr/bin/pidstat', '-hdIruw', '-C', os.path.basename(exec_name)]),
                 'iostat' : lambda sdx: (0, ['/usr/bin/iostat', sdx, '-d', '-x', '-y'])}
    # stat_args = {'pidstat' : lambda exec_name: (WAIT_TIME_FOR_LAUNCHING, ['/usr/bin/pidstat', '-hdIruw', '-C', os.path.basename(exec_name))]
    #              'iostat' : lambda sdx: (0, ['/usr/bin/iostat', 'sdx', '-d', '-x'])}
    def __init__(self, log_root, csv_root, stat_cmd, sampling_int):
        assert stat_cmd in ['iostat', 'pidstat']
        self.log_root = log_root
        self.csv_root = csv_root
        self.stat_logpath = None
        self.pstat = None
        self.fd_log = None
        self.stat_cmd = stat_cmd
        self.sampling_interval = sampling_int
        
        if not os.path.isdir(self.log_root):
            os.makedirs(self.log_root)
        if not os.path.isdir(self.csv_root):
            os.makedirs(self.csv_root)

    def set_sampling_interval(self, interval):
        self.sampling_interval = interval

    def start(self, fn_postfix, cmd_arg):
        ''' fn_postfix - log to file prefix
        for iostat, cmd_arg is sdb, sbc; fr pidstat, target cmd full path '''
        delay, cmd = self.stat_args[self.stat_cmd](cmd_arg)
        cmd.append(str(self.sampling_interval))
        print("INFO start sampling %s @ %d sec, cmd=%s" % (self.stat_cmd, self.sampling_interval, cmd))
        if delay > 0:
            time.sleep(delay)
        fn = "%s_%s.log" % (self.stat_cmd, fn_postfix)
        self.stat_logpath = os.path.join(self.log_root, fn)
        try:
            fd_log = open(self.stat_logpath, 'w')
            #pidstr = pidlist.decode('utf-8').replace('\n', ',')
            self.pstat = Popen(cmd, stdout=fd_log, stderr=fd_log, shell=False)
            self.fd_log = fd_log
        except RuntimeError as rte:
            print("WARN caught RuntimeError %s. No pidstat will be captured" % rte)
            self.stop()
        except Exception as ex:
            print("WARN something is wrong: %s. No pidstat will be captured" %  ex)   #sys.exc_info()
            self.stop()
            
    def stop(self, convert=False):
        if self.fd_log: 
            self.fd_log.close()
            self.fd_log = None
        if self.pstat:
            self.pstat.kill()
            self.pstat = None
        if convert and self.csv_root:
            self.convert_log()

    def convert_log(self, csv_root_path=None):
        csv_root = csv_root_path if csv_root_path else self.csv_root 
        cvs_prefix = os.path.join(csv_root, self.stat_cmd, self.stat_logpath[:-4])
        f_2csv = "make_csv_%s" % self.stat_cmd
        if f_2csv in self.__dir__():
            cvsfn = self.__getattribute__(f_2csv)(self.stat_logpath, cvs_prefix)
            print("INFO %s sampling output to %s" % (self.stat_cmd, cvsfn))
        else:
            print("WARN could not find converter to csv")
        return cvsfn
