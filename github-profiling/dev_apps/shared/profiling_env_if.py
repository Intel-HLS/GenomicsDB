#! /usr/bin/python3

def get_callsets_root():
    raise NotImplementedError("specify the root directory that contains callsets file")

def get_cfg_json_ws():
    raise NotImplementedError("specify the root directory that contains json configuration files")

def get_runlog_root():
    raise NotImplementedError("specify the root directory that is used for logging")

def get_mapping_vid_root():
    raise NotImplementedError("specify the root directory that contains vid map file")

def get_reference_root():
    raise NotImplementedError("specify the root directory that contains reference file")

def get_histogram_root():
    raise NotImplementedError("specify the root directory that contains histogram file")

def get_tdb_ws_root(is_ssd=False):
    raise NotImplementedError("specify the tiledb working directory")
 
def get_interest_pos_root(name=None):
    raise NotImplementedError("specify the root directory for interesting positions files")
##
def get_hosts():
    raise NotImplementedError("specify the host list")
