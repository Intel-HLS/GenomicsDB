---
--- create_scheme.sql
---
--- sqlite3 genomicsdb_loader.db < create_scheme.sql
---- .read pre_fill.sql
---

PRAGMA foreign_keys = ON;
DROP TABLE IF EXISTS loader_config_tag;
DROP TABLE IF EXISTS template;
DROP TABLE IF EXISTS host;

CREATE TABLE IF NOT EXISTS host (
   _id INTEGER PRIMARY KEY AUTOINCREMENT,
   hostname TEXT NOT NULL,
   ipaddress TEXT NULL,
   availibility INTEGER DEFAULT 1,
   creation_ts TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL,
   CONSTRAINT host_uniq UNIQUE (hostname ) ON CONFLICT REPLACE
);
--- template --
--- my_type: 1-vid, 2-callsets, 3-vcf_header, 4- ref_genome, 5=>histogram, 6=>loader, 7=>query
---- params for callsets is default {"ub_callset_row_idx" : 999}
--- TODO add loader_config and query_config
CREATE TABLE IF NOT EXISTS template (
   _id INTEGER PRIMARY KEY AUTOINCREMENT,
   name Text NOT NULL,
   my_type INTEGER DEFAULT 0,
   file_path Path NOT NULL,
   params TEXT NULL,
   extra TEXT NULL,
   CONSTRAINT name_uniq UNIQUE (name ) ON CONFLICT REPLACE
 );
---
--- loader config tags, with default_value
 CREATE TABLE IF NOT EXISTS loader_config_tag (
    _id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    my_type TEXT NOT NULL,
    default_value TEXT NULL,
    tag_code TEXT NULL DEFAULT "",
    CONSTRAINT name_uniq UNIQUE (name ) ON CONFLICT REPLACE
);
--- 
CREATE TABLE IF NOT EXISTS pre-post_task_def (
    _id INTEGER PRIMARY KEY AUTOINCREMENT,
    function_name: TEXT,
    arg_names: TEXT,
    output: TEXT,
    description: TEXT,
    CONSTRAINT name_uniq UNIQUE (function_name) ON CONFLICT REPLACE
);
-- output: success, or errors
-- phases: 1=>pre-load, 2=>post-load, 3=>both
CREATE TABLE IF NOT EXISTS pre_post_loader_task (
    _id INTEGER PRIMARY KEY AUTOINCREMENT,
    task_id REFERENCES pre_post_task_def(_id) NOT NULL,
    batch_id REFERENCES batch_def(_id) NOT NULL,
    inputs: TEXT,
    phases: INTEGER DEFAULT 3,
    description: TEXT,
    CONSTRAINT name_uniq UNIQUE (name) ON CONFLICT REPLACE
);
----
-- define a run, status = run through, canceled, ..?
-- currently only profile with time
-- target_command: cmd_root_path,...
-- loader_configs: for loader, ldnames, for query list of loader_def_id
--- tiledb_ws_root require for loading
--- bound_host: optional, if defined, will run that box
--- input_config: the json let for GenomicsDB loader or query config  
---- 
CREATE TABLE IF NOT EXISTS batch_def (
   _id INTEGER PRIMARY KEY AUTOINCREMENT,
   name TEXT,
   cmd_root TEXT NOT NULL,
   tiledb_ws_root TEXT,
   on_host TEXT NULL,
   mpirun_def TEXT,
   input_config TEXT,
   description TEXT,
   creation_ts TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL,
   CONSTRAINT name_ts_uniq UNIQUE (name) ON CONFLICT REPLACE
);
--- for 'use'' key 
---run_parent: 0=>no, !0=>yes
CREATE TABLE IF NOT EXISTS use_batch_def (
   _id INTEGER PRIMARY KEY AUTOINCREMENT,
   parent_batch_id REFERENCES batch_def(_id),
   run_parent INTEGER 1,
   test_def_list TEXT,
   description TEXT,
   creation_ts TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL
);

----
-- loader_batch_id and query_batch_id mutual exclusive
-- query_config: 
CREATE TABLE IF NOT EXISTS test_def (
  _id INTEGER PRIMARY KEY AUTOINCREMENT,
  batch_id REFERENCES batch_def(_id) NOT NULL,
  num_parallel INTEGER DEFAULT 1,
  lcname TEXT DEFAULT '',
  tiledb_ws TEXT NULL,
  query_config TEXT,
  creation_ts TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL
);
-- default profiling are [time, pidstat]
-- overrite_ws = false is not supported yet
CREATE TABLE IF NOT EXISTS run_def (
  _id INTEGER PRIMARY KEY AUTOINCREMENT,
  hostname TEXT NOT NULL,
  run_defs TEXT NOT NULL,
  overrite_ws INTEGER DEFAULT 1,
  profilers TEXT,
  description TEXT,
  creation_ts TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL
);
--- for analysis, add partition 1 and total size for convenience
CREATE TABLE IF NOT EXISTS test_result (
   _id INTEGER PRIMARY KEY AUTOINCREMENT,
   test_def_id REFERENCES test_def (_id) NOT NULL,
   run_def_id REFERENCES run_def (_id) NOT NULL,
   cmd_version TEXT,
   time_result TEXT NOT NULL,
   target_stdout_path TEXT,
   target_stderr_path TEXT,
   partition_1_size INTEGER DEFAULT 0,
   db_size INTEGER DEFAULT 0,
   pidstat_path TEXT,
   other_profiling TEXT,
   genome_result TEXT,
   comment TEXT,
   creation_ts TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL
);
----
--- execution info, buid_tag is an abitrary value
----
CREATE TABLE IF NOT EXISTS exec_info (
   _id INTEGER PRIMARY KEY AUTOINCREMENT,
   buid_tag TEXT NOT NULL,
   name TEXT NOT NULL,
   version TEXT DEAFULT '0.4',
   full_path TEXT NOT NULL,
   build_time INTEGER NULL,
   hash TEXT,
   other TEXT,
   creation_ts TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL
 );
----
--- operation is 'create' or 'upgrade'
----
CREATE TABLE IF NOT EXISTS scheme_info (
   _id INTEGER PRIMARY KEY AUTOINCREMENT,
   version TEXT DEAFULT '0.2',
   operation TEXT DEFAULT 'create',
   creation_ts TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL
 );
