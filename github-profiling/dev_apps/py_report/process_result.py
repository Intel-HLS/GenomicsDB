#! /usr/bin/python3
#! pylint: disable=broad-except, too-many-instance-attributes
''' u_parse_plot_result2.py, similar to u_parse_plot_result.py
    only work with new group profiling
'''
import os
import os.path
import re
#from time import localtime, strftime
from collections import OrderedDict, namedtuple
#import logging
import glob

import numpy as np
import pandas as pd
from pygenomicsdblib.utils.common import make_csv_iostat
from .result_parser import DensityParser
from .generate_plot import PlotGenerator

DEFAULT_RESULT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../run_logs"))
CommandToken = "Command being timed"
pd.set_option('display.width', 160)
pd.get_option('display.width')

class BaseResultHandler():
    PlotConfig = namedtuple('PlotConfig', ['func', 'x_label', 'y_label'])
    wc_in_sec = lambda st: sum(x * int(t) for x, t in zip([60, 1, 0.001], re.split(r"[:.]", st)))

    utime_keys = [('Elapsed (wall clock) time', lambda ln: ln.split('):')[1].strip(), wc_in_sec, 'Wall Clock in Second'), \
    ("Percent of CPU this job got", lambda ln: ln.split('got:')[1].strip()[:-1], lambda x: int(x)/100, "CPU%")]
    __make_title = lambda self, nn, exp: "Profiling Results of %s at %s%s" % ('_'.join(exp.split('_')[1:]), nn[0].upper(), nn[1:])
    node_name = lambda self, fn: os.path.basename(fn).split('_')[0]

    def __init__(self, proj_name, root_path=DEFAULT_RESULT_ROOT):
        self.root_path = root_path if root_path else DEFAULT_RESULT_ROOT
        self.project_name = proj_name
        self.sel_option = 0             # default plot Elapse time
        self.plotor = PlotGenerator()
        self.parser = DensityParser()
        self.df_stat_dir = None
        self.df_stat_mean = None
        self.f_eid_name = None
        self.f_gid_name = None

    def _get_stats_dirinfo(self, host_list):
        test_root = os.path.join(self.root_path, self.project_name)
        assert os.path.isdir(test_root), '%s not found' % test_root
        samplings = []
        for root, dirs, _ in os.walk(test_root, topdown=True):
            for dname in dirs:
                if dname == "sampling":
                    hostname = os.path.basename(root)
                    if hostname in host_list:
                        test_id = int(os.path.basename(os.path.dirname(root)).split('_')[0][2:])
                        samplings.append([int(test_id/10), test_id%10, os.path.join(root, dname), hostname])
        self.df_stat_dir = pd.DataFrame(samplings, columns=['gid', 'eid', 'sample_root', 'hostname'])

    @staticmethod
    def show_data(name, df):
        for g in df.groupby(['gid']):
            print('%s, 0=', (name, g[0]))     #MING
            print(g[1].head(3))      # MING

    def __prep4plot_data(self, data_set, idx=None, plot_pos_list=None):
        # flat out and convert time to sec
        if plot_pos_list:
            flat_data = [[d, p, s, o, [self.utime_keys[self.sel_option][2](x) for x in v]] for (d, p, s, o), v in data_set.items() if p in plot_pos_list]
        else:
            flat_data = [[d, p, s, o, [self.utime_keys[self.sel_option][2](x) for x in v]] for (d, p, s, o), v in data_set.items()]
        df = pd.DataFrame(flat_data, columns=['density', 'pos', 'seg_size', 'qopt', 'result'])
        if idx is not None:
            df['file_idx'] = idx
        df['wc_mean'] = df.result.apply(np.mean)
        df['wc_std'] = df.result.apply(np.std)
        return df

    def _build_stats_mean_df(self, host_list=None):
        if self.df_stat_dir is None or host_list:
            self._get_stats_dirinfo(host_list)
        mean_list = []
        miss_stat = []
        for _, row in self.df_stat_dir.iterrows():
            os.chdir(row.sample_root)
            for fn in glob.glob("*.log"):
                num_pos = int(fn.split('_')[3])
                src = os.path.join(row.sample_root, fn)
                density = fn.split('_')[2]
                mean_item = {'gid' : row.gid, 'eid' : row.eid, 'density' : density, "num_pos" : num_pos}
                try:
                    f_csv = make_csv_iostat(src, src[:-4])
                    df = pd.DataFrame.from_csv(f_csv).filter(items=["r/s", "rkB/s", "avgrq-sz", "avgqu-sz", "r_await", "%util"])
                    mean_item.update(df.mean().to_dict())
                    mean_list.append(mean_item)
                except Exception:
                    miss_stat.append("%d%d::%s" % (row.gid, row.eid, os.path.basename(src)))
    # old    means = [get_mean_values(f, row) for f in glob.glob("*.log") if int(f.split('_')[3]) in num_pos]
    #        mean_list.extend(means)
        df_stat = pd.DataFrame(mean_list)
# caller reorder df_stat = df_stat[['gid', 'eid', 'density', 'num_pos', "r/s", "rkB/s", "avgrq-sz", "avgqu-sz", "r_await", "%util"]]
        self.df_stat_mean = df_stat
        if miss_stat:
            print("INFO no stat available:", miss_stat)
    def find_results(self):
        '''
        The collected results files *.result
        return: a list contains result file pathes
        '''
        test_root = os.path.join(self.root_path, self.project_name)
        assert os.path.isdir(test_root), '%s not found' % test_root
        results = []

        for root, _, files in os.walk(test_root, topdown=True):
            for fn in files:
                # print('fn=', fn)
                if fn.endswith('result'):
                    results.append(os.path.join(root, fn))
        return results

    def __extract_from_log(self, result_fn):
        ''' one time cmd token at time.
            return a dict: time_cmd_token4exact : [elapse_times | CPU% |.. ]
            parsed_token for "variant-density" is [density, num_pos, seg_size]
        '''
        tm_token = self.utime_keys[self.sel_option][0]
        with open(result_fn, 'r') as ifd:
            # print("INFO: processing ", os.path.basename(result_fn))
            profiled_data = OrderedDict()
            for ln in ifd:
                ln = ln.lstrip()
                if ln.startswith(CommandToken):
                    ln = ln.split(':')[1].strip()
                    key_token = self.parser.parse_command_line(ln[1:-1])  #time_cmd token not include ""
                    if key_token:
                        if key_token not in profiled_data:
                            profiled_data[key_token] = []
                        for lnt in ifd:
                            if lnt.lstrip().startswith(tm_token):
                                time_token = self.utime_keys[self.sel_option][1](lnt)
                                profiled_data[key_token].append(time_token)
                                break
        return profiled_data

    def get_result_df(self, sel_hosts):
        ''' test_name = "variant_density"
            group_list: [10, 20, ..], NSY'''
        assert len(sel_hosts) > 0
        all_result_files = self.find_results()
        df_results = pd.DataFrame([[int(li[0][2]), int(li[0][3]), li[1], fn] for li, fn in [(fn.split(os.sep)[-3:], fn) for fn in all_result_files]], columns=['gid', 'eid', 'hostname', 'fname'])
        df_results['file_idx'] = df_results.index
        # print(df_results.groupby(['hostname']).size())
        df_data = None
        for idx, row in df_results.iterrows():
            profiled_data = self.__extract_from_log(row['fname'])
            df_row = self.__prep4plot_data(profiled_data, idx)
            df_data = df_row if idx == 0 else df_data.append(df_row)
        df_m = pd.merge(df_results, df_data, on='file_idx', how='inner')
        if not isinstance(sel_hosts, list) or sel_hosts not in df_m.hostname.unique():
            print('get_result_df, invalid host name, will select all', sel_hosts)  #MING
            sel_hosts = df_m.hostname.unique()
        df_ret = df_m[(df_m.hostname.isin(sel_hosts))]
        return df_ret


    @staticmethod
    def _get_stat_mean_df(df_samples, num_pos):
        mean_list = []
        for _, row in df_samples.iterrows():
            os.chdir(row.sample_root)
            for fn in glob.glob("*.log"):
                if int(fn.split('_')[3]) in num_pos:
                    src = os.path.join(row.sample_root, fn)
                    density = fn.split('_')[2]
                    mean_item = {'gid' : row.gid, 'eid' : row.eid, 'density' : density}
                    try:
                        f_csv = make_csv_iostat(src, src[:-4])
                        df = pd.DataFrame.from_csv(f_csv).filter(items=["r/s", "rkB/s", "avgrq-sz", "avgqu-sz", "r_await", "%util"])
                        mean_item.update(df.mean().to_dict())
                        mean_list.append(mean_item)
                    except Exception as ede:
                        print("Info exception while process %d%d_%s: %s" % (row.gid, row.eid, os.path.basename(src), ede))
        return mean_list

    def compare_group_by_eids(self, group_eid_pairs, pos_list=None, option=0):
        raise NotImplementedError

    def plot_group_eids(self, eid_list, pos_list=None, option=0):
        raise NotImplementedError
