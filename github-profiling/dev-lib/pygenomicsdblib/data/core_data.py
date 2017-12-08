#! /usr/bin/python3

import os
import os.path
import sqlite3
from utils.common import get_my_logger 

class DataHandler(object):
    CreateDBScript = 'create_scheme.sql'
    PrefillDBScript = 'pre_fill.sql'

    queries = {
        "AllHost" : 'SELECT hostname FROM host;',
        "Host" : 'SELECT hostname FROM host WHERE availibility = 1;',
        "Template" : 'SELECT name, file_path, params, extra FROM template;',
        "LC_Tag" : 'SELECT name, type, default_value FROM loader_config_tag where user_definable=0;',
        "LC_OverrideTag" : 'SELECT name, type, default_value,tag_code FROM loader_config_tag where user_definable=1;',
        'AllUser_LC' : 'SELECT name, config, _id from loader_config_def;',

        'INSERT_LOADER' : "INSERT INTO loader_config_def (name, config) VALUES (\"%s\", \"%s\");",
        'INSERT_BATCH_LOADER_DEF' : 'INSERT INTO batch_def (name, loader_configs, target_cmd_path, tiledb_ws_root, mpirun_def, runtime_env_id) VALUES (\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%s);',
        'INSERT_BATCH_QUERY_DEF' : 'INSERT INTO batch_def (name, loader_configs, target_cmd_path, pick_mode, config_args, runtime_env_id) VALUES (\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%s);',
        'INSERT_TEST_DEF' : 'INSERT INTO test_def (batch_id, num_parallel, full_cmd, tiledb_ws, lcname, query_config) VALUES (%d, %d,\"%s\",\"%s\",\"%s\",\"%s\");',
        'INSERT_RUN_DEF' : 'INSERT INTO run_def (hostname, run_defs) VALUES (\"%s\",\"%s\");',
        
        # batch select
        'ALL_LOADER_BATCH' : 'SELECT bf._id, bf.name, bf.target_cmd_path, bf.loader_configs, bf.tiledb_ws_root, bf.mpirun_def, re.library_path, bf.description FROM batch_def bf LEFT JOIN runtime_env re ON bf.runtime_env_id=re._id  WHERE  tiledb_ws_root is not null;', 
        'ALL_QUERY_BATCH' : 'SELECT bf._id, bf.name, bf.target_cmd_path, bf.loader_configs, bf.tiledb_ws_root, bf.mpirun_def, re.library_path, bf.description FROM batch_def bf LEFT JOIN runtime_env re ON bf.runtime_env_id=re._id  WHERE  tiledb_ws_root is null;',

        'BATCH_4RETESTS' : 'SELECT td._id, td.lcname, td.num_parallel, bd.tiledb_ws_root, td.full_cmd FROM test_def td, batch_def bd WHERE td.batch_id=%d and bd._id = td.batch_id; ',

        'BATCH_TESTS' : 'SELECT td.id FROM test_def td, batch_def bd  WHERE bd.name=\"%s\" and td.batch_id= bd._id; ',

        'User_Config' : 'SELECT config FROM loader_config_def WHERE name in (%s);',
        'User_Config_dict' : "SELECT name, config FROM loader_config_def WHERE name in (%s);",
        'Test_Results' : 'SELECT tr.time_result, tr.genome_result, tr.pidstat_path, rl.lcname, rl.num_parallel,rl.profiler FROM test_result tr, test_def rl where tr.run_id=rl._id and rl.run_def_id=%d order by rl._id desc;',

        'Get_Command' : 'SELECT target_cmd_path from batch_def WHERE _id = %d', 
        #Run_ConfigNames => Batch_COnfigNames
        'Run_ConfigNames' : 'SELECT lcname from test_def WHERE batch_id = %d',
        'Last_Batch_Run_Def' : 'SELECT _id FROM batch_def ORDER BY _id DESC LIMIT 1',
        'Get_LoaderRunId' : 'SELECT loader_configs from batch_def WHERE _id = %d',
        # Target Version
        'SELECT_LAST_TAG' : 'SELECT buid_tag FROM exec_info ORDER BY _id desc limit 1;',
        'SELECT_EXEC_INFO' : 'SELECT name, version, full_path, hash from exec_info where buid_tag=\"%s\";',
        "GET_CMD_VERSION" : "SELECT version, hash FROM exec_info WHERE full_path=\"%s\" ORDER BY _id desc limit 1;",
        "GET_COUNT_BY_HASH" : "SELECT count(*) FROM exec_info WHERE full_path=\"%s\" AND hash=\"%s\";",
        'INSERT_EXEC_INFO' : 'INSERT INTO exec_info (buid_tag, name, version, full_path, build_time, hash) VALUES (\"%s\",\"%s\",\"%s\",\"%s\",%d,\"%s\");',
        
        # env table
        'SELECT_ADDITIONAL_LIBS' : 'SELECT _id, library_path FROM runtime_env;',
        'ADD_ADDITIONAL_LIBS' : 'INSERT INTO runtime_env (library_path) VALUES (\"%s\");',
        
        # update, TODO, not used yet
        'SET_HOSTS' : 'UPDATE host set avalability=%d WHERE hostname in (\"%s\");',

        # remote 
        "GET_MY_TIDS" : "SELECT run_defs FROM run_def WHERE _id=%d;",
        "SELECT_MY_RUNS" : "SELECT td._id, td.full_cmd, td.tiledb_ws, bd.target_cmd_path, re.library_path \
            FROM test_def td, run_def rd, batch_def bd LEFT JOIN runtime_env re ON bd.runtime_env_id=re._id \
            WHERE td._id in (%s) AND rd._id=%d and td.batch_id=bd._id;",
        "INSERT_RESULTS" : "INSERT INTO test_result (test_def_id, run_def_id, cmd_version, time_result, \
            genome_result, partition_1_size, db_size, pidstat_path)  \
            VALUES (%d, %d, \"%s\", \"%s\", \"%s\", %s, \"%s\", \"%s\");" 
    }
    @staticmethod
    def __to_null(x):
        return x if x else 'null'

    def __create_db(self, db_path):
        root_dir = os.path.dirname(__file__)
        creator = os.path.join(root_dir, self.CreateDBScript)
        filler = os.path.join(root_dir, self.PrefillDBScript)
        if os.path.isfile(creator) and os.path.isfile(filler):
            with open(creator) as fdsc:
                sqltext = fdsc.read()
            db_conn = sqlite3.connect(db_path)
            mycursor = db_conn.cursor()
            mycursor.executescript(sqltext)
            with open(filler) as fdsc:
                sqltext = fdsc.read()
            mycursor.executescript(sqltext)
            db_conn.commit()
            self.db_conn = db_conn
        else:
            raise RuntimeError("Miss DB creation scripts")

    def __init__(self, db_path, autoCreateDB=False):
        self.logger = get_my_logger(self.__class__.__name__)
        if os.path.isfile(db_path):
            self.db_path = db_path
            self.db_conn = sqlite3.connect(self.db_path)
        elif autoCreateDB:
            try:
                self.__create_db(db_path)
                self.db_path = db_path
            except:
                raise RuntimeError("cannot create database")
        else:
            raise NameError("Invalid database name")

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __del__(self):
        self.close()

    def get_db_path(self):
        return self.db_path if self.db_path else None

    def close(self):
        if self.db_conn:
            self.db_conn.close()

    def addBatchLoaderConfig(self, name, lcnames, cmd, tdbws, mpirun_def, envid=None):
        mycursor = self.db_conn.cursor()
        val = self.__to_null(envid)
        stmt = self.queries['INSERT_BATCH_LOADER_DEF'] % (name, lcnames, cmd, tdbws, mpirun_def, val)
        mycursor.execute(stmt)
        batch_id = mycursor.lastrowid
        self.db_conn.commit()
        mycursor.close()
        return batch_id

    def addBatchQueryConfig(self, name, loader_cfgs, cmd, pick_mode, mycfgs, envid=None):
        mycursor = self.db_conn.cursor()
        val = self.__to_null(envid)
        stmt = self.queries['INSERT_BATCH_QUERY_DEF'] % (name, loader_cfgs, cmd, pick_mode, mycfgs, val)
        mycursor.execute(stmt)
        batch_id = mycursor.lastrowid
        self.db_conn.commit()
        mycursor.close()
        return batch_id

    def getBatch4ReTests(self, batch_id):
        ''' get tests list of a batch '''
        assert batch_id, 'batch_id is required'
        mycursor = self.db_conn.cursor()
        query = self.queries["BATCH_4RETESTS"] % batch_id
        ret = []
        for row in mycursor.execute(query).fetchall():
            ret.append((row[0], row[1], row[2], row[3], row[4]))
        mycursor.close()
        return ret

    def getBatchTests(self, batch_name):
        assert batch_name, 'batch_name is required'
        mycursor = self.db_conn.cursor()
        query = self.queries["BATCH_TESTS"] % batch_name
        return [row[0] for row in mycursor.execute(query)]

    def getAllLoaderBatches(self):
        mycursor = self.db_conn.cursor()
        for row in mycursor.execute(self.queries['ALL_LOADER_BATCH']):
            yield row[0], row[1], row[2], row[3], row[6], row[7], [row[4], row[5]]
        mycursor.close()
        
    def getAllQueryBatches(self):
        mycursor = self.db_conn.cursor()
        for row in mycursor.execute(self.queries['ALL_QUERY_BATCH']):
            yield row[0], row[1], row[2], row[3], row[6], row[7], [row[4], row[5]]
        mycursor.close()

    def getTheseRunsInfo(self, runid, loader_list):
        assert runid, 'run id is required'
        mycursor = self.db_conn.cursor()
        query = self.queries["Runs_of_RunDef"] % runid
        run_info = []
        for row in mycursor.execute(query):
            if row[0] in loader_list:
                run_info.append(dict({'lc_name': row[0], 'num_proc': row[1], 'tdb_ws':row[2], 'host':row[3], 'loader_config':row[4].split()[-1]}))
        return run_info

    def getRunConfigs(self, runid, bFillFlag=True):
        ''' fillFlag = 0 => no fill, = 1 => fill with '-'; = 2 => fill with default_value '''
        mycursor = self.db_conn.cursor()
        definable_tags = self.getUserDefinableConfigTags(mycursor)
        myrunid = runid if runid else self.getLastRun()
        cfg_names = self.getRunConfigNames2(myrunid, mycursor)
        if cfg_names:
            query = self.queries['User_Config'] % ",".join(["\"%s\"" % x for x in cfg_names])
            lc_items = []
            for row in mycursor.execute(query):
                lc_items.append(dict(row[0]))     # user overrided tags
            
            if bFillFlag > 0:                           # need to fill not overrided tags
                for cfg in lc_items:
                    for key, val in definable_tags.items():
                        if key not in cfg:
                            cfg[key] = str(val[2])
                        else:
                            cfg[key] = "%s*" % cfg[key]
            return myrunid, lc_items
        else:
            return None, None

    def getLastRun(self):
        mycursor = self.db_conn.cursor()
        last_run_def_id = mycursor.execute(self.queries['Last_Batch_Run_Def']).fetchone()
        mycursor.close()
        return last_run_def_id[0]

    def getRunCommand(self, runid):
        mycursor = self.db_conn.cursor()
        row = mycursor.execute(self.queries['Get_Command'] % runid).fetchone()
        cmd = row[0] if row else None
        mycursor.close()
        return cmd

    # for review ??
    def getRunConfigNames2(self, runid, cursor=None):
        ''' return the loader config list for a run, last run if runid is None  '''
        mycursor = cursor if cursor else self.db_conn.cursor()
        ret = []
        stmt = self.queries['Run_ConfigNames'] % runid
        for row in mycursor.execute(stmt):
            ret.append(row[0]) 
        if not cursor:
            mycursor.close()
        return ret

    def getLoadRunId(self, batchid):
        mycursor = self.db_conn.cursor()
        runid = mycursor.execute(self.queries['Get_LoaderRunId'] % batchid).fetchone()
        return list(runid[0])

    def getRunConfigsDict(self, myrunid, bFillFlag=True):
        ''' fillFlag = 0 => no fill, = 1 => fill with '-'; = 2 => fill with default_value '''
        assert myrunid, 'run id is required'
        mycursor = self.db_conn.cursor()
        definable_tags = self.getUserDefinableConfigTags(mycursor)

        cfg_names = self.getRunConfigNames2(myrunid, mycursor)
        if cfg_names:
            query = self.queries['User_Config_dict'] % ",".join(["\"%s\"" % x for x in cfg_names])
            lc_items = {}
            for row in mycursor.execute(query):
                lc_items[row[0]] = dict(row[1])     # user overrided tags
            
            if bFillFlag > 0:                # need to fill not overrided tags
                for cfg in lc_items.values():
                    for key, val in definable_tags.items():
                        if key not in cfg:
                            cfg[key] = str(val[2])
                        else:
                            cfg[key] = "%s*" % cfg[key]
            return lc_items
        else:
            return None
    def getAllResult(self, runidstr):
        assert runidstr
        mycursor = self.db_conn.cursor()
        all_results = []
        runid = int(runidstr)
        for row in mycursor.execute(self.queries['Time_Results'] % runid):
            try:
                rowresult = dict()
                rowresult['rtime'] = dict(row[0])
                rowresult['gtime'] = row[1]
                rowresult['pidstat'] = row[2]
                rowresult['lcname'] = row[3]
                rowresult['n_parallel'] = row[4]
                rowresult['extra_info'] = dict(row[5]) if row[5] else None
            except:
                print("getAllResult() exception: ignored item for runid", runidstr)

            all_results.append(rowresult)
        cmd = mycursor.execute(self.queries['Get_Command'] % runid).fetchone()
        return all_results, cmd[0]

    def getHosts(self):
        mycursor = self.db_conn.cursor()
        mycursor.execute(self.queries["Host"])
        rows = list(mycursor.fetchall())
        hosts = []
        for r in rows:
            hosts.append(r[0])
        mycursor.close()
        return hosts

    def getTemplates(self):
        mycursor = self.db_conn.cursor()
        for r in mycursor.execute(self.queries['Template']):
            yield r[0], r[1], r[2], r[3] 
        mycursor.close()

    def getUserDefinableConfigTags(self, cursor):
        mycursor = cursor if cursor else self.db_conn.cursor()
        lc_overridable_tags = {}
        for row in mycursor.execute(self.queries['LC_OverrideTag']):
            lc_overridable_tags[row[0]] = list(row)
        if not cursor:
            mycursor.close()
        return lc_overridable_tags

    def getConfigTags(self):
        lc_fixed_tags = {}
        mycursor = self.db_conn.cursor()
        for row in mycursor.execute(self.queries['LC_Tag']):
            lc_fixed_tags[row[0]] = list(row)
        lc_overridable_tags = self.getUserDefinableConfigTags(mycursor)
        mycursor.close()
        return lc_fixed_tags, lc_overridable_tags         
    
    def getAllUserDefinedConfigItems(self):
        defined_loaders = {}
        mycursor = self.db_conn.cursor()
        for row in mycursor.execute(self.queries['AllUser_LC']):
            line = row[1].replace("u'", "\"").replace("'", "\"")
            cfg = eval(line)
            defined_loaders[row[0]] = (cfg, row[2])  
        mycursor.close()
        return defined_loaders
    
    def addUserDefinedConfig(self, lcname, configStr):
        mycursor = self.db_conn.cursor()
        stmt = self.queries['INSERT_LOADER'] % (lcname, configStr)
        mycursor.execute(stmt)
        loader_id = mycursor.lastrowid
        self.db_conn.commit()
        mycursor.close()
        return loader_id

    def addBatchTest(self, batch_id, cmd, tiledb_ws, loader_cfg, num_parallel=1, query_cfg=None):
        mycursor = self.db_conn.cursor()
        stmt = self.queries['INSERT_TEST_DEF'] % (batch_id, num_parallel, cmd, tiledb_ws, loader_cfg, query_cfg)
        mycursor.execute(stmt)
        test_id = mycursor.lastrowid
        self.db_conn.commit()
        mycursor.close()
        return test_id  

    # target version
    def hasVersionInfo(self, full_fn, hashcode):
        stmt = self.queries['GET_COUNT_BY_HASH'] % (full_fn, hashcode)
        mycursor = self.db_conn.cursor()
        ret = False
        ret = mycursor.execute(stmt).fetchone()[0] > 0
        mycursor.close()
        return ret           
    def addVersionInfo(self, tag, target_name, version, full_fn, build_time, hash_code):
        stmt = self.queries["INSERT_EXEC_INFO"] % (tag, target_name, version, full_fn, build_time, hash_code)
        #print("get_exec_info: stmt=", stmt)
        mycursor = self.db_conn.cursor()
        mycursor.execute(stmt)
        run_id = mycursor.lastrowid
        self.db_conn.commit()
        mycursor.close()
        return run_id 

    def checkExecsGen(self):
        mycursor = self.db_conn.cursor()
        my_tag = mycursor.execute(self.queries['SELECT_LAST_TAG']).fetchone()
        if my_tag:
            stmt = self.queries['SELECT_EXEC_INFO'] % my_tag
            for row in mycursor.execute(stmt):
                yield row[1], row[2], row[3]
            mycursor.close()
        else:
            mycursor.close()
            raise RuntimeError("No version info in db")
    
    def get_last_version_info(self, full_path):
        mycursor = self.db_conn.cursor()
        stmt = self.queries['GET_CMD_VERSION'] % full_path
        row = mycursor.execute(stmt).fetchone()
        mycursor.close()
        return row if row else None

    # additional env info
    def getAdditionalLibs(self):
        ret = {}
        mycursor = self.db_conn.cursor()
        stmt = self.queries['SELECT_ADDITIONAL_LIBS']
        for row in mycursor.execute(stmt):
            ret[row[0]] = row[1]
        mycursor.close()
        return ret

    def addAdditionalLibs(self, libpath):
        mycursor = self.db_conn.cursor()
        stmt = self.queries['ADD_ADDITIONAL_LIBS'] % libpath
        mycursor.execute(stmt)
        envid = mycursor.lastrowid
        self.db_conn.commit()
        mycursor.close()
        return envid

    def addRunDef(self, hostname, test_list):
        mycursor = self.db_conn.cursor()
        stmt = self.queries['INSERT_RUN_DEF'] % (hostname, test_list)
        mycursor.execute(stmt)
        rid = mycursor.lastrowid
        self.db_conn.commit()
        mycursor.close()
        return rid

    def getMyRuns(self, task_id):
        mycursor = self.db_conn.cursor()
        stmt = self.queries['GET_MY_TIDS'] % task_id
        tids = mycursor.execute(stmt).fetchall()[0][0]
        ret = []
        stmt = self.queries['SELECT_MY_RUNS'] % (tids, task_id)
        #  print("INFO %s: getMyRuns stmt=%s" %(task_id, stmt))
        for row in mycursor.execute(stmt):
            ret.append( (row[1], row[2], row[0], row[3], row[4]) )
        mycursor.close()
        return ret

    def saveResults(self, task_id, test_id, exec_info, time_results, genome_results, cvsfiles, tdb_size):
        mycursor = self.db_conn.cursor()
        stmt = self.queries['INSERT_RESULTS'] % (task_id, test_id, exec_info, str(time_results), str(genome_results), 0, tdb_size, cvsfiles)
        print(stmt)
        mycursor.execute(stmt)
        rid = mycursor.lastrowid
        self.db_conn.commit()
        mycursor.close()
        return rid
        