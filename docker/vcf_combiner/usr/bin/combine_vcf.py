#! /usr/bin/python3
# pylint: disable=missing-docstring, invalid-name, broad-except, too-many-branches, too-many-locals, line-too-long

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
from subprocess import Popen
import json
from getopt import getopt
import logging
import uuid
from collections import namedtuple
from collections import OrderedDict
import types
import vcf

Version = '0.6'
ProductName = os.path.basename(__file__)

COL_PARTITION_SIZE_UNIT = 16384
DefaultVIDFile = "/usr/share/cont-intel/vid.json"
logging.basicConfig(stream=sys.stdout, level=logging.INFO)

##### loader config file tags
loader_cfg_t = {
    "produce_combined_vcf": True,
    "num_parallel_vcf_files": 1,
    "callset_mapping_file" : "",
    "vcf_output_format" : "z",
    "vcf_output_filename" : "",
    "index_output_VCF": True,
    "reference_genome" : "",
    "column_partitions" : [],
    "do_ping_pong_buffering" : True,
    "offload_vcf_output_processing" : True,
    "produce_GT_field" : False,
    "size_per_column_partition" : COL_PARTITION_SIZE_UNIT,
    "vid_mapping_file": ''
}
cp_pos = namedtuple('pos_only', "begin, vcf_output_filename")
cp_chr = namedtuple('chromosome', 'begin, end, vcf_output_filename')
##### loader config file tags

def get_loader_cfg(**kwargs):
    ''' named tuple loader config '''
    f = lambda x: x() if isinstance(x, types.FunctionType) else x
    updated = {k : f(kwargs[k]) if k in kwargs else v for k, v in loader_cfg_t.items()}
    np_LoaderCfg = namedtuple('loader_cfg', ','.join(updated.keys()))
    return np_LoaderCfg(**updated)

def get_col_partition(output_fn, begin, chromosome=None, end=None):
    cp_nt = cp_pos(begin=begin, vcf_output_filename=output_fn) if not chromosome else \
        cp_chr(begin={chromosome:begin}, end={chromosome:end}, vcf_output_filename=output_fn)
    return cp_nt._asdict()

class CombineVCFException(Exception):
    pass

class CombineVCF(object):
    ''' VCF file combiner '''
    logger = logging.getLogger("CombineVCF")
    brief_options = "i:o:R:c:p"
    full_options = ['samples=', 'output=', 'reference=', 'callsets=', 'vid_mapping_file=', \
        'produce_GT_field', 'chromosome=', 'begin=', 'end=', 'dryrun']

    def __init__(self):
        self.dryrun = False
        self.output_file = None
        self.vid_mapping_file = DefaultVIDFile

    def _parse_args(self, args):
        def check_chromosome():
            assert begin, 'No begin position is given'
            if end:
                assert end >= begin, 'End position must be greater or equal to begin position'
            return get_col_partition(self.output_file, begin, chromosome, end if end else None)

        myopts, _ = getopt(args, self.brief_options, self.full_options)
        callset_mapping_file = None
        produce_GT_field = False
        chromosome = None
        begin = None
        end = None
        vid_mapping_file = DefaultVIDFile
        for opt, user_input in myopts:
            if opt == '-p' or opt == '--produce_GT_field':
                produce_GT_field = True
            elif opt == '--dryrun':
                self.dryrun = True
            else:
                assert user_input, 'specify a value for option %s' % opt
                if opt == '-i' or opt == '--samples':
                    file_list = user_input.split(',')
                    vcf_inputfiles = self.__get_inputs(file_list)
                elif opt == '-o' or opt == '--output':
                    self.output_file = self.__check_output(user_input)
                elif opt == '-R' or opt == '--reference':
                    assert os.path.isfile(user_input) or os.path.islink(user_input), "specify a valid reference file name"
                    reference_genome = user_input
                elif opt == '-c' or opt == '--callsets':
                    assert os.path.isfile(user_input), "specify a valid callset file name"
                    callset_mapping_file = user_input
                    num_part_units = self.__check_callset(callset_mapping_file)
                elif opt == '--vid_mapping_file':
                    assert os.path.isfile(user_input), "specify a valid vid mapping file"
                    vid_mapping_file = user_input
                elif opt == '--chromosome':
                    chromosome = user_input
                elif opt == '--begin':
                    begin = int(user_input)
                elif opt == '--end':
                    end = int(user_input)
                else:
                    print("WARN: unknown option %s, ignored", opt)

        if not callset_mapping_file:
            callset_mapping_file = "callsets_%s.json" % datetime.now().strftime("%y%m%d%H%M")
            num_part_units = self.__generate_callsets_json(vcf_inputfiles, callset_mapping_file)

        assert reference_genome, "missing reference file"
        assert self.output_file, "missing output file"
        assert num_part_units != 0, "No valid callset/sample found in input files"

        col_par_setting = check_chromosome() if chromosome \
            else get_col_partition(self.output_file, begin if begin else 0)
        loader_cfg = get_loader_cfg()
        loader_cfg = loader_cfg._replace(reference_genome=reference_genome, \
            vcf_output_filename=self.output_file, \
            column_partitions=[col_par_setting],  \
            size_per_column_partition=int(abs(num_part_units)) * COL_PARTITION_SIZE_UNIT, \
            callset_mapping_file=callset_mapping_file, \
            vid_mapping_file=vid_mapping_file, \
            produce_GT_field=True if produce_GT_field else False)
        return loader_cfg

    @staticmethod
    def __check_output(user_input):
        if not os.path.exists(os.path.dirname(user_input)):
            os.makedirs(os.path.dirname(user_input))
        else:
            assert not os.path.isdir(user_input), "specify a valid output file name"
        return user_input

    @staticmethod
    def __check_callset(callset_mapping_file):
        cs_cfg = json.load(callset_mapping_file)
        return 0 - len(cs_cfg["callsets"])

    def __generate_callsets_json(self, vcf_inputfiles, json_fname):
        global_callset_idx = 0
        callsets_dict = OrderedDict()
        for vcf_file in vcf_inputfiles:
            read_mode = 'r' if vcf_file[-3:] == 'vcf' else 'rb'
            with open(vcf_file, read_mode) as fd:
                vcf_reader = vcf.Reader(fd)
                local_callset_idx = 0
                for callset_name in vcf_reader.samples:
                    curr_callset_info = OrderedDict()
                    if callset_name in callsets_dict:
                        ss = str(uuid.uuid4())
                        self.logger.warning('Duplicate callset name %s: appending _%s', callset_name, ss)
                        callset_name += ('_' + ss)
                    curr_callset_info["row_idx"] = global_callset_idx
                    curr_callset_info["idx_in_file"] = local_callset_idx
                    curr_callset_info["filename"] = vcf_file
                    callsets_dict[callset_name] = curr_callset_info
                    local_callset_idx += 1
                    global_callset_idx += 1
        with open(json_fname, 'w') as ofd:
            json.dump({'callsets' : callsets_dict}, ofd, indent=4, separators=(',', ': '))
        return global_callset_idx

    @staticmethod
    def __is_vcf_file_list(fn):
        if fn[-3:] == '.gz' or fn[-4:] == 'vcf':
            return False
        with open(fn, 'r') as fd:
            for line in fd:
                line = line.strip()
                if line:
                    if fn[-3:] == '.gz' or fn[-4:] == 'vcf' or os.path.isfile(line):
                        return True
        return False

    def __get_inputs(self, inputfiles):
        if len(inputfiles) == 1 and self.__is_vcf_file_list(inputfiles[0]):
            with open(inputfiles[0], 'r') as fd:
                inputs = [line.strip() for line in fd if os.path.isfile(line.strip())]
        else:
            inputs = [line.strip() for line in inputfiles if os.path.isfile(line)]
        return inputs if inputs else RuntimeError("No valid samples input files found")

    def generate_loader_config(self, nt_loader):
        json_fname = os.path.join(os.path.dirname(self.output_file), \
            "loader_config_%s.json" % datetime.now().strftime("%y%m%d%H%M"))
        with open(json_fname, 'w') as ofd:
            json.dump(nt_loader._asdict(), ofd)
        return json_fname

    @staticmethod
    def combine(loader_config):
        the_exec_cmd = ['vcf2tiledb', loader_config]
        pexec = Popen(the_exec_cmd, shell=False)
        pexec.communicate()
        return pexec.wait()

    def run(self):
        err = None
        try:
            nt_loader_cfg = self._parse_args(sys.argv[1:])
            loader_cfg = self.generate_loader_config(nt_loader_cfg)
            if self.dryrun:
                self.logger.info("Done dryrun: generated loader_cfg at " + loader_cfg)
            else:
                ret_code = self.combine(loader_cfg)
                if ret_code != 0:
                    raise CombineVCFException("ERROR: failed to combine vcf files. Error code = " + str(ret_code))
                elif os.path.isfile(self.output_file):
                    self.logger.info("Done - combined vcf sample files into a single vcf file " + self.output_file)
                else:
                    raise CombineVCFException("Internal error. Could not produce combined VCF file from vcf sample files")
        except CombineVCFException as myex:
            self.logger.exception(myex)
            raise
        except Exception as ex:
            self.logger.exception(ex)
            err = sys.exc_info()[1]
        finally:
            if err:
                raise CombineVCFException("Failed to combining VCF files: %s" % (err))

    @staticmethod
    def get_my_name():
        return ProductName

def test_code_runs_after_pylint():
    ''' a quick test for github check-in '''
    print("Run a quick test ...")
    x = CombineVCF().get_my_name()
    assert x
    assert x == ProductName

if __name__ == "__main__":
    combiner = CombineVCF()
    if len(sys.argv) > 1:
        combiner.run()
    else:
        combiner.get_my_name()
