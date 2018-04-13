# /usr/bin/python3
'''
python3 -m unittest utest_arg -c -f     # no .py
python3 utest_arg.py -c -f              #

-c for catching ^C
-f for stop on 1st failure

'''
import sys
import os
import re
import logging
import unittest.mock as mock

logger = logging.getLogger(__name__)

# patch 1) histogram generating,
# 2)  start_pos = HistogramManager(self.params['H'].the_val).calc_bin_begin_pos(self.params['P'].the_val, vb, ve), mock.patch()
# 3) args2 = {'R' : 'path2ref', 'H' : 'path2histogram', 'range' : 'chr1:34:23888445, chr4::32323888445', 'C' : 'path2callsets'}
#TEST_DATA_ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
TEST_DATA_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

PART_START = [1, 222222, 333333, 4444444, 5555555, 66666666, 77777777, 888888888, 999999999, 1000000000]
_callsets_fn = os.path.join(TEST_DATA_ROOT, 'test_resource', 'callsets.json')

test_args = {
    'R' : os.path.join(TEST_DATA_ROOT, 'test_resource', "Homo_sapiens_assembly19.fasta"),
    'C' : [_callsets_fn if os.path.isfile(_callsets_fn) else os.path.dirname(_callsets_fn), os.path.join(TEST_DATA_ROOT, 'test_output')],
    'H' : "/path/to/histogram",
    'i' : os.path.join(TEST_DATA_ROOT, 'test_resource', "vcfs"),
    'o' : os.path.join(TEST_DATA_ROOT, 'test_output'),
    'V' : [os.path.join(TEST_DATA_ROOT, 'test_resource', 'vid.json')],  # "GRCH37" causes fail
    'range' : ['4:1-191154276,1:1-249250621,2:1-243199373', '1-249250621,249250622-492449994', '1:1-249250621,2:1-243199372', '4,6:1-171115066'],
    'P' : len(PART_START),
    'positions' : os.path.join(TEST_DATA_ROOT, 'test_resource', 'query_examples', 'query_1_10_f1.json'),
    'A': ["DP", "MQ"],
    'print-AC': None,
    'print-calls': None
}

arg_names = lambda l: [re.sub(r'\W+', '', x) for x in l]
def print_bad_usage(text, neg_char='-'):
    print(neg_char * 40)
    print(text)
    print()

def print_good_usage(text, pos_char='+'):
    print(pos_char * 40)
    print(text)
    print()

def my_options(my_opt_list):
    return {x : (test_args[x][0] if isinstance(test_args[x], list) else test_args[x]) for x in arg_names(my_opt_list) if x in test_args}

def get_arg(key, idx=0):
    if key in test_args:
        if isinstance(test_args[key], list):
            return test_args[key][idx]
        return test_args[key]
    return None

def get_query_args():
    return my_options(['R', 'C', 'o', 'V', 'positions', 'print-AC'])

def get_mock_env():
    predefined_vid = os.path.join(TEST_DATA_ROOT, 'docker_src', 'known_vid_mapper_files')
    mock_env = mock.Mock()
    mock_env.get_mapping_vid_root.return_value = predefined_vid
    return mock_env

def get_config_fpath(fn):
    return os.path.join(TEST_DATA_ROOT, 'test_output', fn)

def get_ldcfg_4query():
    fname = os.path.join(TEST_DATA_ROOT, 'test_resource', 'test_ld_cfg.json')
    assert os.path.isfile(fname)
    return fname

def get_inputs4callsets():
    callsets_fn = os.path.join(TEST_DATA_ROOT, 'test_output', 'gen_callsets.json')
    logger.debug('TEST_DATA_ROOT=%s", callsets_fn=%s, input=%s', repr(TEST_DATA_ROOT), repr(callsets_fn), test_args['i'])
    return {'C': callsets_fn, 'i' : test_args['i']}

def get_number_part(idx):
    return len(test_args['range'][idx].split(','))

# ALL_ARG_CFG = {"R" : VAL_R, "C" : VAL_C, "V" : VAL_V, "i" : VAL_i, "o" : VAL_o, "H" : VAL_H, "P" : VAL_P, "--range" : VAL_range}
# def arg_tuple(**kwargs):
#     argd = {re.sub(r'\W+', '', k).lower() : v for k,v in kwargs.items()}
#     return namedtuple('test_arg', tuple(argd))(**argd)
#aa = arg_tuple(**args2)

# VAL_R = os.path.join(TEST_DATA_ROOT, 'test_resource', "Homo_sapiens_assembly19.fasta")
# VAL_C = [_callsets_fn if os.path.isfile(_callsets_fn) else os.path.dirname(_callsets_fn), os.path.join(TEST_DATA_ROOT, 'test_output')]
# VAL_H = "/path/to/histogram"
# VAL_i = os.path.join(TEST_DATA_ROOT, 'test_resource', "vcfs")
# VAL_o = os.path.join(TEST_DATA_ROOT, 'test_output')
# VAL_V = [os.path.join(TEST_DATA_ROOT, 'test_resource', 'vid.json')]  # "GRCH37" causes fail
# VAL_range = ['1-3000, 32-50000', '1:1-200000,CHR1:300000-400000', 'CHR5:1-4000000:', 'CH4,CH11:-345444444']
# VAL_P = len(PART_START)
def my_options_old(my_opt_list=None):
    ''' obsolete '''
    # test_dp = importlib.import_module(m_all_options) if m_all_options else sys.modules[__name__]
    test_dp = sys.modules[__name__]
    arg_n = [x for x in dir(test_dp) if x.startswith('VAL_')]
    if my_opt_list:
        arg_n = [x for x in arg_n if x.split('_')[-1] in arg_names(my_opt_list)]
    arg_val = [test_dp.__dict__.get(x) for x in arg_n]
    return dict(zip([x.split('_')[-1] for x in arg_n], arg_val))
