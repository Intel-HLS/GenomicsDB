#! /usr/bin/python3

import json
import os
import os.path
import platform
import shutil
from enum import Enum
#from utils.common import *
from utils.constant_def import GDB_COMMANDS, ENV_TILEDB_ROOT, DEFAULT_TDB_PREFIX

def check_files(lib_files):
    ret = True
    if isinstance(lib_files, list):
        for libfile in lib_files:
            if not os.path.isfile(libfile):
                print(" Not find required tlibrary file", libfile)
                ret = False
    elif not os.path.isfile(lib_files):
        print(" Not find required tlibrary file", lib_files)
        ret = False
    return ret

class Tags(Enum):
    tag_name = "name"
    tag_loader_configs = "loader_configs"
    tag_additional_libs = "additional_libs"
    tag_tiledb_ws_root = "tiledb_ws_root"
    tag_pararell = "pararell_def"
    tag_loader_cmd = "command_path"

class LoaderConfig(object):
    def __init__(self, loader_config):
        with open(loader_config, 'r') as ldf:
            json_def = json.load(ldf)

        assert Tags.tag_name.value in json_def, 'name must be defines'
        assert Tags.tag_loader_cmd.value in json_def, '%s must be defined' % Tags.tag_loader_cmd.value
        self.__check_cmd(json_def[Tags.tag_loader_cmd.value])
        assert Tags.tag_loader_configs.value in json_def, '%s must be defined' \
         % Tags.tag_loader_configs.value
        self.__check_tdb_ws(json_def[Tags.tag_tiledb_ws_root.value] \
         if Tags.tag_tiledb_ws_root.value in json_def else None)
        if Tags.tag_additional_libs.value in json_def:
            if not check_files(json_def[Tags.tag_additional_libs.value]):
                print(" Not all required library files found")

        if Tags.tag_pararell.value in json_def:
            num_config = len(json_def[Tags.tag_loader_configs.value])
            for pos, _ in json_def[Tags.tag_pararell.value].items():
                assert (int(pos) < num_config), 'invalid %s value %d' % (Tags.tag_pararell.value, pos)
        self.json_def = json_def

    
    def __check_tdb_ws(self, tdb_root):
        if not tdb_root:
            tdb_root = os.environ.get(ENV_TILEDB_ROOT, os.path.dirname(DEFAULT_TDB_PREFIX))
        self.tdb_ws_prefix = os.path.join(tdb_root, os.path.basename(DEFAULT_TDB_PREFIX)) 

    def __check_cmd(self, cmd_path):
        cmd_name = os.path.basename(cmd_path)
        assert (cmd_name in GDB_COMMANDS), 'unknown command %s' % cmd_name
        if not os.path.dirname(cmd_path):
            cmd_path = shutil.which(cmd_path)
        if platform.system() != 'Windows':
            assert os.path.exists(cmd_path), "not find required execution file %s " % cmd_path
        self.command_path = cmd_path

    def get_loaders(self):
        return self.json_def[Tags.tag_loader_configs.value]

    def get_name(self):
        return self.json_def[Tags.tag_name.value]

    def get_tileDB_prefix(self):
        return self.tdb_ws_prefix

    def get_additional_libs(self):
        return self.json_def[Tags.tag_additional_libs.value]if Tags.tag_additional_libs.value in self.json_def else None

    def get_mpirun_def(self):
        return self.json_def[Tags.tag_pararell.value] if Tags.tag_pararell.value in self.json_def else None
    
    def get_loader_command(self):
        return self.command_path
