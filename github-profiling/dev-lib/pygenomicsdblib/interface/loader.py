#! /bin/python3
import os
import os.path
from utils.common import get_my_logger, is_windows, get_default_workdir
from utils.constant_def import DEFAULT_DB_NAME, DEFAULT_WORKSPACE
from core import loader_profiler

class BaseInterface(object):

    def __init__(self, working_dir=None):
        self.working_dir = working_dir if working_dir else get_default_workdir()
        self.loader = None
    def get_hostlist(self):
        pass

    def enable_hosts(self, hostlist):
        pass

    def disable_hosts(self, hostlist):
        pass

class LoaderInterface(BaseInterface):
    '''TODO: get loader_config_list(), copy_loader_config between dbs; '''
    loaders = {}

    def __init__(self, working_dir=None):
        BaseInterface.__init__(self, working_dir)
        self.logger = get_my_logger(self.__class__.__name__)
        self.is_dryrun = is_windows()
        dbpath = os.path.join(self.working_dir, DEFAULT_DB_NAME)
        # if os.path.exists(dbpath):
        self.use_database(dbpath)

    def use_database(self, db_path):
        assert db_path, 'path to database cannot be empty'
        prf_key = os.path.abspath(db_path)
        if prf_key not in self.loaders:
            self.loaders[prf_key] = loader_profiler.LoaderProfiler(db_path)
        self.loader = self.loaders[prf_key]
        
    def __check_config_file(self, config_file):  #was __parse__config
        ''' parse loader config file '''
        assert self.loader, "use_database first"
        if not os.path.isfile(config_file):
            config_file = os.path.join(self.working_dir, config_file)
            assert os.path.isfile(config_file),  'invalid config file %s' % config_file
        return config_file

    def run(self, config_file, host_list=None, dryrun=None):
        ''' from config file '''
        assert self.loader, "use_database first"
        config_file = self.__check_config_file(config_file)
        if dryrun or self.is_dryrun:
            self.loader.dry_run(config_file, host_list) 
        else:
            self.loader.launch(config_file, host_list) 
    
    def generate_batch_loader_config(self, config_file, root_path=None):
        assert self.loader, "use_database first"
        if not root_path:
            root_path = os.path.join(self.working_dir, config_file[:-4])
        files = self.loader.generate_batch_loader_config(config_file, root_path)
        print("Generated @ %s: %s" % (root_path, files))

    def run_batch(self, batch_name, host_list=None, dryrun=None):
        ''' batch by name '''
        assert self.loader, "use_database first"
        if dryrun or self.is_dryrun:
            self.loader.dry_run_batch(batch_name, host_list) 
        else:
            self.loader.launch_batch(batch_name, host_list) 

    def run_tests(self, test_id_list, host_list=None, dryrun=False):
        ''' adhoc run,  hostname=None => use assigned '''
        assert self.loader, "use_database first"
        if dryrun:
            self.loader.dry_run_tests(test_id_list, host_list) 
        else:
            self.loader.launch_tests(test_id_list, host_list)

    def close_all(self):
        for ldr in self.loaders.values():
            ldr.close()
        self.loaders.clear()
        if self.loader:
            self.loader = None
