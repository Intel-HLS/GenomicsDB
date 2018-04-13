#! /usr/bin/python3

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
from subprocess import CalledProcessError
import logging
from collections import namedtuple
import re
import json
from .metadata_generator import generate_callsets_json, generate_histogram

logger = logging.getLogger(__name__)

ParamHandler = namedtuple('ParamHandler', ['the_val', 'default_val', 'handler', "is_generate", 'help_msg'])

is_vcf_file = lambda fn: (fn[-3:] == '.gz' or fn[-4:] == 'vcf') and fn[0] != '.'

class DockerParamHandlerBase:
    params = None
    known_vid_mapper = {}
    LoadConfigCheckTags = [("column_partitions", None), ("callset_mapping_file", 'C'), ("vid_mapping_file", 'V'), ("reference_genome", 'R')]

    def check_load_config(self, apath):
        with open(apath, 'r') as ifd:
            cs_data = json.load(ifd)
        for x in self.LoadConfigCheckTags:
            assert x[0] in cs_data, "Invalid load config file. miss %s" % x[0]
        return cs_data

    @staticmethod
    def _generate_hist(to_file, aparams):
        assert 'V' in aparams and 'C' in aparams, "missing vid mapping file or callsets"
        try:
            generate_histogram('vcf_histogram', aparams['V'].the_val, aparams['C'].the_val, to_file)
        except CalledProcessError:
            raise RuntimeError("cannot find vcf_histogram")

    @staticmethod
    def _generate_callsets(to_file, aparams):
        assert 'i' in aparams, "missing input VCFs"
        _ = generate_callsets_json(aparams['i'].the_val, to_file)

    def _check_gen_path(self, apath, mykey):
        logger.debug('check_gen_path, param=', apath)
        if os.path.isdir(apath):
            gen_file_path = os.path.join(apath, self.params[mykey].default_val)
            self.params[mykey].handler(gen_file_path, self.params)
            self.params[mykey] = self.params[mykey]._replace(the_val=gen_file_path, is_generate=True)
        else:
            self.params[mykey] = self.params[mykey]._replace(the_val=apath)

    def _check_inputs(self, user_input, mykey):
        inputfiles = user_input.split(',')
        vcf_files = None
        if len(inputfiles) == 1:
            apath = inputfiles[0]
            if os.path.isdir(apath):
                vcf_files = [os.path.abspath(os.path.join(apath, fn)) for fn in os.listdir(apath) if is_vcf_file(fn)]
            else:
                with open(apath, 'r') as fd:
                    vcf_files = [line.strip() for line in fd if os.path.isfile(line.strip())]
        else:
            vcf_files = [os.path.abspath(line.strip()) for line in inputfiles if os.path.isfile(line.strip())]

        assert vcf_files, "not valid input VCF files found. Usage: " + self.params[mykey].help_msg
        self.params[mykey] = self.params[mykey]._replace(the_val=vcf_files)

    def _load_predefined_vid(self):
        vid_root = self._myenv.get_mapping_vid_root()
        if os.path.isdir(vid_root):
            return {os.path.splitext(f)[0] : f for f in os.listdir(vid_root) if f[-4:] == '.vid'}        # print("predefined vid root is ", vid_root)
        print("invalid path to predefined vid files %s" % vid_root)
        return None

    def _check_vid(self, apath, mykey):
        if os.path.isfile(apath):
            self.params[mykey] = self.params[mykey]._replace(the_val=os.path.abspath(apath))
        else:
            if not self.known_vid_mapper and self._myenv:
                self.known_vid_mapper = self._load_predefined_vid()
            if self.known_vid_mapper and apath in self.known_vid_mapper:
                fn_path = os.path.join(self._myenv.get_mapping_vid_root(), self.known_vid_mapper[apath])
                assert os.path.isfile(fn_path), "cannot find predefined vid mapper file %s, please contact GenomicsDB team" % apath
                self.params[mykey] = self.params[mykey]._replace(the_val=fn_path)
            else:
                raise RuntimeError("Unknown vid mapper file " + apath)

    def _check_output_dir(self, apath, mykey):
        # if os.path.exists(apath):
        #     shutil.rmtree(apath)
        # os.makedirs(apath)
        open(os.path.join(apath, '__tiledb_workspace.tdb'), 'w').close()
        self.params[mykey] = self.params[mykey]._replace(the_val=apath)

    def __init__(self, option_inputs, myenv):
        self._myenv = myenv
        self.params = dict()
        option_list = [re.sub(r'\W+', '', p) for p in option_inputs]

        if 'R' in option_list:
            self.params['R'] = ParamHandler("", None, None, False, 'the path to a reference_genome file')
        if 'C' in option_list:
            self.params['C'] = ParamHandler("", 'callsets.json', self._generate_callsets, False, 'the path to a callsets file or a root path where generated callsets will be stored')
        if 'H' in option_list:
            self.params['H'] = ParamHandler("", 'histogram', self._generate_hist, False, 'the path to a histogram file or a root path where generated histogram will be stored. Required for vcf_importer')
        if 'V' in option_list:
            self.params['V'] = ParamHandler("", None, None, False, 'the path to a vid mapping file or a name of predefined vid mapping files')
        if 'i' in option_list:
            self.params['i'] = ParamHandler("", None, None, False, 'root path of input VCF files or path to a VCF file list or the option followed by a set of VCF files separated by ","')
        if 'o' in option_list:
            self.params['o'] = ParamHandler("", None, None, False, 'root path of GenomicsDB workspace')
        if 'P' in option_list:
            self.params['P'] = ParamHandler(-1, 1, None, False, "number of partitions, defulat is 1")

    def validate(self, my_args):
        raise NotImplementedError("abstract method.")

    def get_meta_files(self):
        assert self.params and 'R' in self.params and 'V' in self.params and 'C' in self.params
        return self.params['R'].the_val, self.params['V'].the_val, self.params['C'].the_val
