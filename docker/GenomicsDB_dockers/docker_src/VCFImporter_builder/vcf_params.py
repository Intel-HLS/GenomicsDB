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
import os
import os.path
import json
import re
import logging
from pygenomicsdbperf.utils.cfg_json_templates import LOADER_CFG_T, COL_PARTITION_SIZE_UNIT
from pygenomicsdbperf.shared.histogram import HistogramManager
from pygenomicsdbperf.shared.docker_param_base import DockerParamHandlerBase, ParamHandler

logger = logging.getLogger(__name__)

PartKeys = [["array", "begin", "end", "workspace"], ["array", "begin", "workspace"]]
# PartMaker = lambda n, c, b, e, w: (PartKeys[0], ["TEST%d" % n, dict([(c, b)]) if b else c, dict([(c, e)]) if e else c, w]) if c else ((PartKeys[0], ["TEST%d" % n, b, e, w]) if e else (PartKeys[1], ["TEST%d" % n, b, w]))
PartMaker = lambda n, c, b, e, w: (PartKeys[0], ["TEST%d" % n, dict([(c, b)]) if b else c, dict([(c, e)]) if e else c, w]) if c else ((PartKeys[0], ["TEST%d" % n, b, e, w]) if e else (PartKeys[1], ["TEST%d" % n, b, w]))

class VCFDockerParamHandler(DockerParamHandlerBase):
    def _no_chromosome_name(self, vb, ve):
        assert self.params and 'H' in self.params and 'P' in self.params
        self._check_gen_path(self.params['H'].the_val, 'H')
        ws = self.get_tiledb_ws()
        start_pos = HistogramManager(self.params['H'].the_val).calc_bin_begin_pos(self.params['P'].the_val, vb, ve)
        return [{"array": "TEST%d" % (i), "begin": begin, "workspace": ws} for i, begin in enumerate(start_pos)]

    def _check_partition(self, args, mykey):
        def parse_fields(fields):
            assert fields
            cname = None
            if len(fields) == 1:
                pos_range = fields[0].split('-')
                if len(pos_range) == 1:          # 'chr1'
                    return pos_range[0], 1, None
            else:
                cname = fields[0]
                pos_range = fields[1].split('-')
            vbegin = int(pos_range[0]) if pos_range[0] else 1
            vend = int(pos_range[1]) if pos_range[1] else None
            if vbegin and vend:
                assert vbegin < vend
            return cname, vbegin, vend
        partitions = []
        myws = self.get_tiledb_ws()
        for i, arg in enumerate(args.split(',')):
            cn, vb, ve = parse_fields(arg.split(':'))
            # if i == 0 and not cn:
            #     self._no_chromosome_name(vb, ve)
            #     break
            # elif not cn:
            #     raise ValueError("missing chromosome name")
            # if cn in arg_list:
            #     raise ValueError("duplicated chromosome names")
            partitions.append(dict(zip(*PartMaker(i, cn, vb, ve, myws))))
        self.params[mykey] = self.params[mykey]._replace(the_val=partitions)

    def __init__(self, option_inputs, myenv):
        assert option_inputs, 'empty option list is not allowed'
        option_list = [re.sub(r'\W+', '', p) for p in option_inputs]
        DockerParamHandlerBase.__init__(self, option_list, myenv)
        # ?? super().__init__(option_list, myenv)
        if "range" in option_list:
            self.params['range'] = ParamHandler("", "::", self._check_partition, False, "the RANGE is in format 'chromesome_name:begin:end', separated by ',', all fields are optional. The defaults is full geno. Fields: 1) 'chromesome_name' is a name of chromosome, default is all chromesome; 2) 'begin' is the begining position, default is 1; 'end' is the last position, default is the last position within the chromosome or all for empty chromesome_name. Example: 1) '--range chr1:20:12345678' 2) '--range chr1:20:12345678,chr3,chr6::5678123'")
        # print('Params: ', self.params, option_list)

    def validate(self, my_args):
        if 'R' in self.params:
            assert os.path.isfile(my_args.R), self.params['R'].help_msg
            self.params['R'] = self.params['R']._replace(the_val=os.path.abspath(my_args.R))
        if 'V' in self.params:
            self._check_vid(my_args.V, 'V')
        if 'i' in self.params and my_args.i:
            self._check_inputs(my_args.i, 'i')
        if 'o' in self.params and my_args.o:
            self._check_output_dir(os.path.abspath(my_args.o), 'o')
        if 'C' in self.params:
            self._check_gen_path(os.path.abspath(my_args.C), 'C')
        if "range" in self.params and my_args.range:
            self.params['range'].handler(my_args.range, 'range')

    def get_tiledb_ws(self):
        assert self.params and 'o' in self.params
        return self.params['o'].the_val

    def get_num_partitions(self):
        assert self.params
        if 'range' in self.params:
            return len(self.params['range'].the_val)
        elif 'P' in self.params:
            return self.params['P'].the_val
        else:
            return None

    def get_num_inputs(self):
        assert self.params and 'C' in self.params
        cs_fn = self.params['C'].the_val
        with open(cs_fn, 'r') as ifd:
            cs_data = json.load(ifd)
        if "file_division" in cs_data:
            num_file = len(cs_data["file_division"][0])
        else:
            print("cannot find 'file_division' in callset file ", cs_fn, ", will set number file as 10")
            num_file = 10
        # print('get_num_input, num file=', num_file)
        return num_file

    def create_config(self, to_file):
        ref_geno, vid_mapping, callsets = self.get_meta_files()
        lb, ub = 0, self.get_num_inputs()

        nt_loader = LOADER_CFG_T()
        nt_loader_cfg = nt_loader._replace(column_partitions=self.params['range'].the_val, \
        callset_mapping_file=callsets, vid_mapping_file=vid_mapping, \
        reference_genome=ref_geno, delete_and_create_tiledb_array=False, \
        lb_callset_row_idx=lb, ub_callset_row_idx=ub-1, size_per_column_partition=(ub-lb)*COL_PARTITION_SIZE_UNIT)
        with open(to_file, 'w') as ofd:
            json.dump(nt_loader_cfg._asdict(), ofd, sort_keys=True, indent=4)
