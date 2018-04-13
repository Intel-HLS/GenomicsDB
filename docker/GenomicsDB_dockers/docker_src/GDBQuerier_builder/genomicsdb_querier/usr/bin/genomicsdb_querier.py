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
import shutil
import argparse
import logging
import traceback

from pygenomicsdbperf.interfaces.profiling_env_if import ProfilingEnv
from pygenomicsdbperf.shared.parsing_utils import parse_inputs
from pygenomicsdbperf.utils.common import is_docker
from query_params import QueryParamChecker

Version = '0.5'
ProductName = os.path.splitext(os.path.basename(__file__))[0]
OUTPUT_CONFIG_FN = "docker_querycfg_%s.json"
ArgDescription = "Query GenomicsDB"
MY_VOLUME_PATH = "/myvolume"

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger("ProductName")

required_args = ["-R", "-C", "-V", "-o", "--positions"]
optional_args = ["-A"]
passthrough_args = ['--print-AC', '--print-calls']

if __name__ == "__main__":
    myformat = '%(levelname)s %(module)s:%(lineno)d: %(message)s'
    logging.basicConfig(stream=sys.stdout, format=myformat)
    logging.getLogger().setLevel(os.environ.get('LOGLEVEL', logging.INFO))
logger = logging.getLogger(ProductName)

my_parent = lambda: os.path.dirname(__file__) if os.path.dirname(__file__) else os.getcwd()

def exec_cmd(config_fpath, extra_cmd=""):
    cmd_line = ['gt_mpi_gather', '-j', config_fpath, '-s', str(10000)]
    if extra_cmd:
        cmd_line += extra_cmd
    logger.info('Run: %s', ' '.join(cmd_line))
    return subprocess.call(cmd_line, shell=False)

def get_product_name():
    return "%s %s" % (ProductName, Version)

def get_option_list():
    return required_args + optional_args

def save2volume(config_fpath):
    fopath = os.path.join(MY_VOLUME_PATH, os.path.basename(config_fpath))
    shutil.copyfile(config_fpath, fopath)
    logger.debug("query config file path is @ %s", fopath)

class MyEnv(ProfilingEnv):
    def get_cfg_json_ws(self, cust_name=None):
        thepath = os.path.join("/", "tmp", "test_output") if is_docker() else os.path.join(os.path.dirname(os.path.dirname(my_parent())), "test_output")
        if not os.path.exists(thepath):
            os.makedirs(thepath)
        return thepath

def pre_parse(theparser):
    theparser.add_argument('--print-AC', type=bool, nargs='?', const=True, default=True)
    theparser.add_argument('--print-calls', type=bool, nargs='?', const=True, default=False)

def get_passthrough_cmd(myargs):
    more_cmd = []
    if myargs.print_AC:
        more_cmd.append('--print-AC')
    if myargs.print_calls:
        more_cmd.append('--print-calls')
    return more_cmd

def process(cfg_fpath, myenv=None):
    myparser = argparse.ArgumentParser(description=ArgDescription, prefix_chars='-')
    arg_checher = QueryParamChecker(required_args + optional_args)
    myargs = parse_inputs(myparser, arg_checher, required_args, optional_args, pre_parse)
    arg_checher.validate(myargs)
    more_cmd = get_passthrough_cmd(myargs)

    # dryrun, more_cmd = parse_inputs(myparser, arg_checher, required_args, optional_args, PassthroughHandler)
    if cfg_fpath:
        arg_checher.create_config(cfg_fpath)
        if os.path.exists(MY_VOLUME_PATH):
            save2volume(cfg_fpath)
    return myargs.dry_run, more_cmd

if __name__ == "__main__":
    try:
        config_fn = os.path.join(MyEnv().get_cfg_json_ws(), OUTPUT_CONFIG_FN % datetime.now().strftime("%y%m%d%H%M"))
        dryrun, pt_cmd = process(config_fn)
        if dryrun:
            logger.info("Done, dryrun generated load config file %s", config_fn)
            retcode = 0
        else:
            retcode = exec_cmd(config_fn, pt_cmd)
            logging.debug("Done, loaded VCF files into GenomicsDB return code is %s, the config file at %s", retcode, config_fn)
        exit(retcode)
    except Exception as ex:
        print('Error, caught exception ', ex)
        traceback.print_exc(file=sys.stdout)
        exit(1)
