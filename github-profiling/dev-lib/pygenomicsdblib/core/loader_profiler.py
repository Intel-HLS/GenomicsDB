#! /usr/bin/python3

import os
import os.path
import json
from core.profiler_base import ProfilerBase
from core.loader_config import LoaderConfig
from core.batch_run_mgr import BatchManager
from core.loader_config_mgr import LoaderConfigMgr
from utils.constant_def import MPIRUN_PATH

def loader_decorate(func):
    def wrapper(obj, json_dir, lc_id, tdb_prefix, num_pr=1):
        return func(obj, json_dir, lc_id, tdb_prefix, num_pr)
    return wrapper

class LoaderProfiler(ProfilerBase):
    def __init__(self, data_path, working_dir):
        super(self.__class__, self).__init__(data_path, working_dir)
        self.loader_cfg_manager = LoaderConfigMgr(self.data_handler, self.working_dir)
        self.batch_manager = BatchManager(self.data_handler, 'loader')
        self.tdb_ws = None
        self.target_cmd = None
    
    def __prepare_test(self, ld_def_list):
        '''ld_def_list - list of loader def files or a single loader def file '''
        run_list = []
        loadercfg = LoaderConfig(loader_cfg_json)
        loadercfg.has_preload
        if isinstance(ld_def_list, list):
            for jsondef in ld_def_list:
                runs = self.__prepare_batch(jsondef)
                if runs:
                    run_list.append(runs)
        else:
            run_list = self.__prepare_batch(ld_def_list)
        return run_list
 
    def __prepare_batch(self, loader_cfg_json):
        lcdef_list = self.loader_cfg_manager.add_user_configs(loadercfg.get_loaders())
        if lcdef_list:
            libs = loadercfg.get_additional_libs()
            mpirun_def = loadercfg.get_mpirun_def()
            if mpirun_def:
                mpiruns = {}
                for pos, pnlist in mpirun_def.items():
                    mpiruns[lcdef_list[int(pos)]] = pnlist
            tdb_prefix = loadercfg.get_tileDB_prefix()
            self.target_cmd = loadercfg.get_loader_command()
            name = loadercfg.get_name()
            batch_id, is_new = self.batch_manager.add_config(name, lcdef_list, self.target_cmd, libs, **{'0': tdb_prefix, '1': mpiruns})
            #regenerate *.json from test_def if not new batch
            test_list = self.__regenerate_json(batch_id, name) if not is_new else None
            if not test_list:
                test_list = self.__process_batch(batch_id, name, lcdef_list, tdb_prefix, mpiruns)
            return test_list
        else:
            print("no loader config item found")

    def __regenerate_json(self, batch_id, name):
        json_dir = self.get_test_json_input_path('loader', batch_id, name)
        results = self.data_handler.getBatch4ReTests(batch_id)
        test_list = []              #test id and cmd
        if results:
            for tid, lcname, num_pr, tdb_ws, _ in results:
                json_fn = self.__generate_json(json_dir, lcname, tdb_ws, num_pr) 
                print(' -- generated json file %s' %json_fn)
                test_list.append(tid)
        return test_list

    @loader_decorate
    def __generate_json(self, json_dir, lc_id, tdb_prefix, num_pr):
        self.tdb_ws = "%s_%d-%s/" % (tdb_prefix, num_pr, lc_id)
        load_config, _ = self.loader_cfg_manager.gen_load_config(lc_id, self.tdb_ws) 
        jsonfn = os.path.join(json_dir, "%d-%s.json" % (num_pr, lc_id))
        with open(jsonfn, 'w') as ofd:
            json.dump(load_config, ofd)
        return jsonfn
    
    @loader_decorate
    def __generate_test_def(self, json_dir, lc_id, tdb_prefix, num_pr):
        jsonfn = self.__generate_json(json_dir, lc_id, tdb_prefix, num_pr)
        full_cmd = "%s -np %d %s %s" % (MPIRUN_PATH, num_pr, self.target_cmd, jsonfn) \
        if num_pr > 1 else "%s %s" % (self.target_cmd, jsonfn)
        return full_cmd

    def __process_batch(self, batch_id, batch_name, lcdef_list, tdb_prefix, user_mpirun):
        json_dir = self.get_test_json_input_path('loader', batch_id, batch_name)
        my_generator = self.__iterate_batch(json_dir, lcdef_list, tdb_prefix, user_mpirun, self.__generate_test_def)
        test_list = []              #test id and cmd
        for lc_id, num_pr, full_cmd in my_generator:
            testid = self.data_handler.addBatchTest(batch_id, full_cmd, self.tdb_ws, lc_id, num_pr)
            test_list.append(testid)
        return test_list

    def __iterate_batch(self, json_dir, lcdef_list, tdb_prefix, user_mpirun, handler):
        for lc_id in lcdef_list:
            if lc_id in self.loader_cfg_manager.get_user_loader_cfg():
                if user_mpirun and lc_id in user_mpirun:
                    for num_pr in user_mpirun[lc_id]:
                        ret_data = handler(json_dir, lc_id, tdb_prefix, num_pr)
                        yield lc_id, num_pr, ret_data 
                else:
                    ret_data = handler(json_dir, lc_id, tdb_prefix)
                    yield lc_id, 1, ret_data

    def generate_batch_loader_config(self, loader_cfg_json, root_path=None):
        ret = {}
        loadercfg = LoaderConfig(loader_cfg_json)
        lcdef_list = self.loader_cfg_manager.add_user_configs(loadercfg.get_loaders())
        if lcdef_list:
            mpirun_def = loadercfg.get_mpirun_def()
            if mpirun_def:
                mpiruns = {}
                for pos, pnlist in mpirun_def.items():
                    mpiruns[lcdef_list[int(pos)]] = pnlist
            tdb_prefix = loadercfg.get_tileDB_prefix()
            if not root_path:
                root_path = self.get_test_json_input_path('loader', 0, loadercfg.get_name())
            my_generator = self.__iterate_batch(root_path, lcdef_list, tdb_prefix, mpiruns, self.__generate_json)
            for lc_id, _, jsonfn in my_generator:
                ret[lc_id] = jsonfn
        return ret

    def dry_run(self, jsondef, host_list=None):
        ''' dry run, not call run_exec '''
        run_list = self.__prepare_test(jsondef)
        self.start_dry_run(run_list, host_list)

    def dry_run_batch(self, batch_name, host_list=None):
        '''run by batch name'''
        run_list = self.data_handler.getBatchTests(batch_name)
        self.start_dry_run(run_list, host_list)

    def dry_run_tests(self, host_list=None, *args):
        assert (len(args) > 0), 'need at least one test id'
        test_id_list = [int(arg) for arg in args]
        self.start_dry_run(test_id_list, host_list)

    def launch(self, jsondef, host_list=None):
        test_id_list = self.__prepare_test(jsondef)
        self.start_launch(test_id_list, host_list)

    def launch_batch(self, batch_name, host_list=None):
        test_id_list = self.data_handler.getBatchTests(batch_name)
        self.start_launch(test_id_list, host_list)

    def launch_tests(self, host_list=None, *args):
        assert (len(args) > 0), 'need at least one test id'
        test_id_list = [int(arg) for arg in args]
        self.start_launch(test_id_list, host_list)
