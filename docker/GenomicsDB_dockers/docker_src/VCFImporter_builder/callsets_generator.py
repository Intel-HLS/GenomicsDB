#! /bin/python3
# pylint: disable=broad-except

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
import re
import argparse
import traceback
import logging
from pygenomicsdbperf.shared.parsing_utils import parse_inputs
from pygenomicsdbperf.shared.docker_param_base import DockerParamHandlerBase
from pygenomicsdbperf.shared.metadata_generator import generate_callsets_json

if __name__ == "__main__":
    myformat = '%(levelname)s %(module)s:%(lineno)d: %(message)s'
    logging.basicConfig(stream=sys.stdout, format=myformat)
    logging.getLogger().setLevel(os.environ.get('LOGLEVEL', logging.DEBUG))
logger = logging.getLogger(__name__)

ArgDescription = "Generates callsets from input VCF files"

required_args = ["-C", "-i"]
optional_args = []

class MyParamHandler(DockerParamHandlerBase):
    def __init__(self, option_inputs, myenv=None):
        assert option_inputs, 'empty option list is not allowed'
        option_list = [re.sub(r'\W+', '', p) for p in option_inputs]
        DockerParamHandlerBase.__init__(self, option_list, myenv)

    def validate(self, my_args):
        if 'i' in self.params and my_args.i:
            self._check_inputs(my_args.i, 'i')
        logger.debug('input=%s, #input=%d, callsets=%s', repr(my_args.i), len(self.params['i'].the_val), repr(my_args.C))
        _ = generate_callsets_json(self.params['i'].the_val, my_args.C)

def process(cfg_fpath=None, theenv=None):
    myparser = argparse.ArgumentParser(description=ArgDescription, prefix_chars='-')
    arg_checher = MyParamHandler(required_args)
    myargs = parse_inputs(myparser, arg_checher, required_args, None)
    arg_checher.validate(myargs)

if __name__ == "__main__":
    try:
        process()
        exit(0)
    except Exception as ex:
        logger.debug('Error, caught exception %s', ex)
        traceback.print_exc(file=sys.stdout)
        exit(1)
