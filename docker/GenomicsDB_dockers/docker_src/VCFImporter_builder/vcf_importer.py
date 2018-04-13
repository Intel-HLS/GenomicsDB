#! /bin/python3
# pylint: disable=missing-docstring, broad-except, abstract-method, too-few-public-methods

"""
  The MIT License (MIT)
  Copyright (c) 2016-2017 Intel Corporation

  Permission is hereby granted, free of charge, to any person obtaining a copy of
  this software and associated documentation files (the "Software"), to deal in
  the Software without restriction, including without limitation the rights to
  use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
  the Software, and to permit persons to whom the Software is furnished to do so,
  subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
  FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
  COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
  IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import sys
import os
import os.path
from datetime import datetime
import subprocess
import argparse
import logging
import shutil
import traceback
from pygenomicsdbperf.interfaces.profiling_env_if import ProfilingEnv
from pygenomicsdbperf.utils.common import is_docker
from pygenomicsdbperf.shared.parsing_utils import parse_inputs

from vcf_params import VCFDockerParamHandler

Version = '0.7'
# the file name reflects the poduct, so one day the skeleton of this file could be generated.
ProductName = os.path.splitext(os.path.basename(__file__))[0]

OUTPUT_CONFIG_FN = "docker_importcfg_%s.json"
ArgDescription = "Process input arguments for %s" % ProductName

# if this volume is mounted, the generated config files will copied there.
MY_VOLUME_PATH = "/myvolume"

# set the required and optional command arguments
required_args = ["-R", "-C", "-V", "-o", "-i"]
optional_args = ["--range"]

if __name__ == "__main__":
    myformat = '%(levelname)s %(module)s:%(lineno)d: %(message)s'
    logging.basicConfig(stream=sys.stdout, format=myformat)
    logging.getLogger().setLevel(os.environ.get('LOGLEVEL', logging.INFO))
logger = logging.getLogger(ProductName)

my_parent = lambda: os.path.dirname(__file__) if os.path.dirname(__file__) else os.getcwd()

class MyEnv(ProfilingEnv):
    def get_mapping_vid_root(self):
        return os.path.join("/", 'usr', 'share', 'cont-intel') if is_docker() else os.path.join(os.path.dirname(my_parent()), "known_vid_mapper_files")

    def get_cfg_json_ws(self, cust_name=None):
        thepath = os.path.join("/", "tmp", "test_output") if is_docker() else os.path.join(os.path.dirname(os.path.dirname(my_parent())), "test_output")
        if not os.path.exists(thepath):
            os.makedirs(thepath)
        return thepath

def get_product_name():
    return "%s %s" % (ProductName, Version)

def get_option_list():
    return required_args + optional_args

def save2volume(config_fpath):
    fopath = os.path.join(MY_VOLUME_PATH, os.path.basename(config_fpath))
    shutil.copyfile(config_fpath, fopath)
    logger.debug("load config file path is @ %s", fopath)

# call the geneomicsDB util command.
def exec_cmd(config_fpath, num_parallel=1):
    cmd_line = ['vcf2tiledb', config_fpath] if num_parallel == 1 else ['/usr/lib64/mpich/bin/mpirun', '-np', str(num_parallel), 'vcf2tiledb', config_fpath]
    logger.info('Run: %s', ' '.join(cmd_line))
    return subprocess.call(cmd_line, shell=False)

def process(cfg_fpath, theenv=None):
    myparser = argparse.ArgumentParser(description=ArgDescription, prefix_chars='-')
    arg_checher = VCFDockerParamHandler(required_args + optional_args, theenv)
    myargs = parse_inputs(myparser, arg_checher, required_args, optional_args)
    arg_checher.validate(myargs)
    if cfg_fpath:
        arg_checher.create_config(cfg_fpath)
        if os.path.exists(MY_VOLUME_PATH):
            save2volume(cfg_fpath)
    return myargs.dry_run, arg_checher.get_num_partitions()

if __name__ == "__main__":
    myenv = MyEnv()
    try:
        config_fn = os.path.join(myenv.get_cfg_json_ws(), OUTPUT_CONFIG_FN % datetime.now().strftime("%y%m%d%H%M"))
        dryrun, num_part = process(config_fn, myenv)
        if dryrun:
            logger.info("Done, dryrun generated load config file %s", config_fn)
            retcode = 0
        else:
            retcode = exec_cmd(config_fn, num_part)
            logging.debug("Done, loaded VCF files into GenomicsDB return code is %s, the config file at %s", retcode, config_fn)
        exit(retcode)
    except Exception as ex:
        print('Error, caught exception ', ex)
        traceback.print_exc(file=sys.stdout)
        exit(1)
