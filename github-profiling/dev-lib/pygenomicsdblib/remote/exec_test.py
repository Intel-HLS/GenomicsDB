#! /usr/bin/python3

import sys
import os
import os.path
import platform
import traceback
from collections import deque
from subprocess import Popen, PIPE, check_output
from data.core_data import DataHandler
from utils.constant_def import TIME_FORMATTER, genome_profile_tags, genome_queryprof_tags, LOADER_EXEC
from remote.get_exec_info import GenomicsExecInfo
from remote.prerun_check import check_tdb_ws, check_libraries, check_executions
from remote.stat_capturer import StatCapturer

DevNull = open(os.devnull, 'wb', 0)

class RunTest():
    def __init__(self, db_path):
        assert os.path.isfile(db_path), 'invalid db path %s' % db_path
        self.db_path = db_path
        self.hostname = platform.node().split('.')[0]

    def __fetch_time_result(self, qexec):
        timeline = qexec.pop().decode('utf-8').strip()
        return dict(x.split('~') for x in timeline.split(','))

    def __fetch_genome_result(self, qexec, line_parser):
        genome_result = []
        for i, gl in enumerate(qexec):
            line = gl.decode('utf-8').strip()
            if line:
                gr = line_parser(line)
                if gr:
                    genome_result.append(gr)
            else:
                print("INFO @%s: empty line %i " % (self.hostname, i))
        return genome_result

    @staticmethod
    def __genome_loader_parser(geno_str):
        lines = geno_str.split(',')
        if lines[0].strip().upper() == 'GENOMICSDB_TIMER':
            ret = {}
            op_str = lines[1].strip().lower()

            if op_str in genome_profile_tags:
                ret['op'] = genome_profile_tags[op_str]
            else:
                ret['op'] = op_str.replace(' ', '_')
            ret['0'] = lines[3]
            ret['1'] = lines[5]
            ret['2'] = lines[7]
            ret['3'] = lines[9]
            ret['4'] = lines[11]
            return ret
        else:
            print("INFO: ignore no GENOMICSDB_TIMER, %s..." % (geno_str[:160]))

    @staticmethod
    def __genome_query_parser(geno_str):
        lines = geno_str.split(',')
        if lines[0].strip().upper() == 'GENOMICSDB_TIMER':
            ret = {}
            op_str = lines[1].strip().lower()

            if op_str in genome_queryprof_tags:
                ret['op'] = genome_queryprof_tags[op_str]
            else:
                print('WARN: operation string %s not found' % (op_str))
                ret['op'] = op_str.replace(' ', '_')
            
            for i in range(2, len(lines)):
                if 'Cpu time(s)' in lines[i]:
                    i += 1
                    ret['0'] = lines[i]    # CPU
                if 'Wall-clock time(s)' in lines[i]:
                    i += 1
                    ret['1'] = lines[i]    # wall clock
            return ret

    @staticmethod
    def get_tdb_size(tiledb_ws):
        return check_output(['du', '-s', tiledb_ws]).decode('utf-8').split()[0]

    @staticmethod
    def launch_time(full_cmd):
        the_exec_cmd = ['/usr/bin/time', "-f", TIME_FORMATTER] + full_cmd
        pexec = Popen(the_exec_cmd, shell=False, stdout=DevNull, stderr=PIPE)
        if not pexec:
            raise RuntimeError("launch test failed")
        return pexec

    def run_test(self, runid, log_path, pid_csv_path):
        #full_cmd,tiledb_ws,_id,target_cmd_path,library_path
        with DataHandler(self.db_path) as dbh:
            task_list = dbh.getMyRuns(runid)
        if not task_list:
            raise RuntimeError("no command found for runid=%d" % (runid))
        stat_capturer = StatCapturer(log_path, pid_csv_path)
        
        for full_cmd, tiledb_ws, tid, target_cmd, libpath in task_list:
            try:
                if libpath:
                    lib_path = libpath if isinstance(libpath, list) else [libpath]
                    check_libraries(lib_path)
                else:
                    check_libraries()
                version_str = GenomicsExecInfo(self.db_path).get_version_info(target_cmd);
                isLoading = LOADER_EXEC == os.path.basename(target_cmd)
                print("++++START @ %s: run_id=%d, test_id=%d, cmd=%s extra_lib=%s; version=%s[%s]" % (self.hostname, runid, tid, target_cmd, libpath, version_str, target_cmd))
                check_executions(target_cmd)
                check_tdb_ws(isLoading, tiledb_ws)
                
                pexec = self.launch_time(full_cmd.split())
                partial_logname = "%d-%d-%s" % (runid, tid, self.hostname)
                stat_capturer.start(partial_logname, target_cmd)
                with pexec.stderr:
                    qexec = deque(iter(pexec.stderr.readline, b''))
                rc = pexec.wait()
                stat_capturer.stop()
                if rc == 0:
                    time_results = self.__fetch_time_result(qexec)
                    line_parser = self.__genome_loader_parser if isLoading else self.__genome_query_parser
                    genome_results = self.__fetch_genome_result(qexec, line_parser)

                    if time_results and genome_results:
                        tdb_size = self.get_tdb_size(tiledb_ws)
                        cvsfiles = stat_capturer.convert_pid_log()
                        with DataHandler(self.db_path) as dbh:
                            dbh.saveResults(runid, tid, version_str, time_results, genome_results, cvsfiles, tdb_size)
                        print("INFO %s: The exit status=%s" % (self.hostname, time_results['9']))
                    else:
                        exit_code = time_results['9'] if time_results else 'unknown'
                        print("WARN %s: No time measured for genomics executable for test %d. The exit status=%s" %  (self.hostname, tid, exit_code))
                else:
                    print("ERROR: error in launching; ", rc)
            except:
                print("  !!!EXCEPT @ %s: run_id=%d, test_id=%d: " % (self.hostname, runid, tid))
                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_exception(exc_type, exc_value, exc_traceback, limit=20, file=sys.stdout)

            finally:
                print("----END @ %s: run_id=%d, test_id=%d," % (self.hostname, runid, tid))
