#! /usr/bin/python3

import os
import os.path
import re
import json
from pygenomicsdbperf.utils.cfg_json_templates import QUERY_CFG_T
from pygenomicsdbperf.shared.docker_param_base import DockerParamHandlerBase, ParamHandler
from pprint import pprint

class QueryParamChecker(DockerParamHandlerBase):
    def _check_positions(self, apath, mykey):
        '''
        position file: [ [ [0, 100 ], 500 ] ]
        '''
        assert os.path.isfile(apath), self.params[mykey].help_msg
        with open(apath, 'r') as ifd:
            pos_list = json.load(ifd)
        assert isinstance(pos_list, list), 'Invalid query position format, see help at https://github.com/Intel-HLS/GenomicsDB/wiki/Querying-GenomicsDB#query-result-format'
        self.params[mykey] = self.params[mykey]._replace(the_val=pos_list)

    # def _from_load_config(self, apath, mykey):
    #     cs_data = self.check_load_config(apath)
    #     # _ = [assert x[0] in cs_data, "Invalid load config file. miss %s" % x[0] for x in lc_tuple]
    #     self.params.update({x[1]: ParamHandler(cs_data[x[0]], None, None, False, None) for x in self.LoadConfigCheckTags if x[1]})
    #     self.params[mykey] = self.params[mykey]._replace(the_val=cs_data[self.LoadConfigCheckTags[0][0]])

    def _check_attributes(self, vals, mykey):
        attr_list = set(self.params[mykey].default_val).union(set(vals))
        self.params[mykey] = self.params[mykey]._replace(the_val=list(attr_list))

    def __init__(self, option_inputs, myenv=None):
        assert option_inputs, 'empty option list is not allowed'
        option_list = [re.sub(r'\W+', '', p) for p in option_inputs]
        DockerParamHandlerBase.__init__(self, option_list, myenv)
        if "positions" in option_list:
            self.params["positions"] = ParamHandler("", None, None, False, 'the path to position list file')
        if 'A' in option_list:
            self.params['A'] = ParamHandler("", ["REF", "ALT", "GT"], None, False, 'the query attributes')

    def validate(self, my_args):
        if 'R' in self.params:
            assert os.path.isfile(my_args.R), self.params['R'].help_msg
            self.params['R'] = self.params['R']._replace(the_val=os.path.abspath(my_args.R))
        if 'V' in self.params:
            self._check_vid(my_args.V, 'V')
        if 'o' in self.params and my_args.o:
            self._check_output_dir(os.path.abspath(my_args.o), 'o')
        if 'C' in self.params:
            self._check_gen_path(os.path.abspath(my_args.C), 'C')
        # if 'L' in self.params:
        #     self._from_load_config(os.path.abspath(my_args.L), 'L')
        if "positions" in self.params:
            self._check_positions(os.path.abspath(my_args.positions), 'positions')
        if 'A' in self.params:
            arg = my_args.A.split(',') if isinstance(my_args.A, str) else my_args.A
            self._check_attributes(arg, 'A')

    def get_tiledb_ws(self):
        assert self.params and 'o' in self.params
        return self.params['o'].the_val
    def get_positions(self):
        assert self.params and 'positions' in self.params
        return self.params['positions'].the_val
    def get_attributes(self):
        assert self.params and 'A' in self.params
        return self.params['A'].the_val

    def create_config(self, to_file):
        ref_geno, vid_mapping, callsets = self.get_meta_files()
        nt_query = QUERY_CFG_T(1)
        nt_cfg = nt_query._replace(callset_mapping_file=callsets, query_attributes=self.get_attributes(), \
        reference_genome=ref_geno, query_column_ranges=self.get_positions(), vid_mapping_file=vid_mapping, \
        workspace=self.get_tiledb_ws())
        with open(to_file, 'w') as ofd:
            json.dump(nt_cfg._asdict(), ofd, sort_keys=True, indent=4)
