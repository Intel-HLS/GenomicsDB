#! /usr/bin/python3
import os
from data.core_data import DataHandler
from utils.common import remake_path, check_x_mode, get_stats_path, get_default_workdir
from utils.constant_def import RUN_SCRIPT, CONFIG_DIRNAME, DEFAULT_WORKSPACE

class ProfilerBase(object):
    def __init__(self, data_path, working_dir):
        self.working_dir = working_dir if working_dir else get_default_workdir()
        self.data_handler = DataHandler(data_path, True)
        self.my_hostlist = self.data_handler.getHosts()

    def verify_host_list(self, host_list):
        ''' verify if host is on our known host_list'''
        ret = [host for host in host_list if host in self.my_hostlist]
        assert ret, 'host %s are not available' % host_list
        return ret

    def close(self):
        self.data_handler.close()

    def get_test_json_input_path(self, test_type, batch_id, name):
        ret = os.path.join(self.working_dir, CONFIG_DIRNAME, "%s-%d-%s" % (test_type, batch_id, name))
        remake_path(ret)
        return ret

    def check_remote_run(self):
        run_script = os.path.join(self.working_dir, RUN_SCRIPT)
        check_x_mode(run_script)
        return run_script

    def __assign_host(self, test_list, romoterun_hosts):
        ''' return {host:run_id} '''
        num_host = len(romoterun_hosts)
        host_runlist = {}
        for i, test_id in enumerate(test_list):
            host = romoterun_hosts[i % num_host]
            if host in host_runlist:
                host_runlist[host].append(test_id) 
            else:
                host_runlist[host] = [test_id]
        run_list = {h : self.data_handler.addRunDef(h, ','.join([str(i) for i in htl])) \
        for h, htl in host_runlist.items()}
        return run_list

    def __prepare_launch(self, test_list, host_list):
        host_list = self.verify_host_list(host_list) if host_list else self.my_hostlist
        host_tasks = self.__assign_host(test_list, host_list) 
        run_script = self.check_remote_run()
        assert run_script, 'cannot find %s at %s ' % (RUN_SCRIPT, self.working_dir)
        statpaths = get_stats_path(self.working_dir)
        for spath in statpaths:
            remake_path(spath)

        return run_script, host_tasks

    def start_dry_run(self, test_list, host_list):
        run_script, host_tasks = self.__prepare_launch(test_list, host_list)
        for host, taskid in host_tasks.items():
            shell_cmd = "ssh %s %s %d %s &" % (host, run_script, taskid, self.data_handler.get_db_path())
            print('DRYRUN: os.system(%s)' % shell_cmd)

    def start_launch(self, test_list, host_list):
        run_script, host_tasks = self.__prepare_launch(test_list, host_list)
        for host, taskid in host_tasks.items():
            lcmd="ssh -tt %s %s %d %s &" % (host, run_script, taskid, self.data_handler.get_db_path())
            print("**launch cmd=", lcmd)
            os.system(lcmd)
            #os.system("ssh -tt %s %s %d %s &" % (host, run_script, taskid, self.data_handler.get_db_path()))
