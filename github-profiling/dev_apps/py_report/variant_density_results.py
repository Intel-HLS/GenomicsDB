#! /usr/bin/python3
#! pylint: disable=broad-except, no-self-use
''' Density Profiling Report
use for auto reload
    %load_ext autoreload
    %autoreload 2

'''
import sys
import traceback
import itertools
from IPython.display import HTML, display
import tabulate
from .process_result import BaseResultHandler

BY_CPU = 1
BY_TIME = 0

def plot_decorator(func):
    def wrapper(obj, *args, **kwargs):
        newkwargs = kwargs if kwargs else dict()
        newkwargs['eid_name'] = obj.f_eid_name
        newkwargs['gid_name'] = obj.f_gid_name
        # print("args=%s plot_deco: newkwargs=%s" % (args, newkwargs))
        return func(obj, *args, **newkwargs)
    return wrapper

class VariantDensity(BaseResultHandler):
    filters = None

    def __init__(self, proj_name, root_path, gidlookup, eidlookup):
        assert eidlookup and gidlookup
        super(VariantDensity, self).__init__(proj_name, root_path)
        self.f_eid_name = lambda x: eidlookup[x] if x < len(eidlookup) and x >= 0 else 'be'+str(x)
        self.f_gid_name = lambda x: gidlookup[x-1] if x <= len(gidlookup) and x > 0 else 'bg'+str(x)
        self.run_name = lambda self, r: "%s and %s mpirun" % (self.f_eid_name(r.eid), self.f_gid_name(r.gid))
        self.sel_host = None

    @staticmethod
    def show_separator():
        display(HTML("<style>hr {border-top: 3px solid navy;}</style><hr>"))

    @plot_decorator
    def f_plot_by_eid(self, df, pos, eid_list, **kwargs):
        ''' run_id: (gid, eid)'''
        # print("f_plot_by_eid: kw=", kwargs)
        assert 'eid_name' in kwargs and 'gid_name' in kwargs
        # eid_list = [rid[1] for rid in rid_list]
        df_list = [df[(df.pos == pos) & (df.eid == eid_list[0])]]
        node = ','.join(df_list[0].hostname.unique())
        if len(eid_list) > 1:
            df_list.append(df[(df.pos == pos) & (df.eid == eid_list[1])])
            title = "Compare %s vs %s for %d positions @ %s" % \
            (kwargs['eid_name'](eid_list[0]), kwargs['eid_name'](eid_list[1]), pos, node)
        else:
            title = "Profiling %s positions for %s @%s" % (pos, kwargs['eid_name'](eid_list[0]), node)
        x_ticklabels = [kwargs['gid_name'](x) for x in df_list[0].gid.unique()]
        return df_list, title, x_ticklabels

    @plot_decorator
    def f_plot3D_by_eid(self, df, pos_list, eid_list, **kwargs):
        ''' run_id: (gid, eid)'''
        # print("f_plot_by_eid: kw=", kwargs)
        assert 'eid_name' in kwargs and 'gid_name' in kwargs
        # eid_list = [rid[1] for rid in rid_list]
        df_list = []
        z_ticklabels = []
        for pos in sorted(pos_list, reverse=True):
            for eid in eid_list:
                df_list.append(df[(df.pos == pos) & (df.eid == eid)])
                z_ticklabels.append("%s, P%d" % (kwargs['eid_name'](eid), pos))
        node = ','.join(df_list[0].hostname.unique())
        title = "Compare %s samples vs %s samples @ %s" % (kwargs['eid_name'](eid_list[0]), kwargs['eid_name'](eid_list[1]), node)
        x_ticklabels = [kwargs['gid_name'](x) for x in df_list[0].gid.unique()]
        return df_list, title, x_ticklabels, z_ticklabels

    @plot_decorator
    def f_plot_by_n_samples(self, df, pos, rid_list, **kwargs):
        ''' run_id: (gid, eid)'''
        # print("f_plot_by_eid: kw=", kwargs)
        assert 'eid_name' in kwargs and 'gid_name' in kwargs
        f_title = lambda x: kwargs['eid_name'](x).split()[0]  # only HDD or SSD
        gids = set([rid[0] for rid in rid_list])
        # print('f_plot_by_n_samples, rid_list=', rid_list, ', gids=', gids)
        assert len(gids) == 1
        gid = gids.pop()
        eid_hdd_list = [rid[1] for rid in rid_list if rid[1] % 2 == 0]
        eid_ssd_list = [rid[1] for rid in rid_list if rid[1] % 2 == 1]
        # print('eid_hdd_list=', eid_hdd_list, ', eid_ssd_list=', eid_ssd_list)
        df_list = []
        if eid_hdd_list and eid_ssd_list:
            assert len(eid_hdd_list) == len(eid_ssd_list)
            df_list.append(df[(df.pos == pos) & (df.gid == gid) & df.eid.isin(eid_hdd_list)])
            df_list.append(df[(df.pos == pos) & (df.gid == gid) & df.eid.isin(eid_ssd_list)])
            title = "%s vs %s for %d pos %s mpirun" % \
            (f_title(eid_hdd_list[0]), f_title(eid_ssd_list[0]), pos, kwargs['gid_name'](gid))
        else:
            id_list = eid_hdd_list if eid_hdd_list else eid_ssd_list
            assert len(id_list) == 1
            df_list.append(id_list)
            title = "%s for %d pos %s mpirun" % (f_title(id_list[0]), pos, kwargs['gid_name'](gid))
        x_ticklabels = [kwargs['eid_name'](x).split()[1] for x in df_list[0].eid.unique()]
        return df_list, title, x_ticklabels

    def compare_group_by_eids(self, group_eid_pairs, pos_list=None, option=0, show_stat=None):
        ''' my_test_name = "variant density"
            group_eid_list = [(0, 1), (2, 3)], high value first. a graph per each (eid1, eid2)
            option = 0 - wall clock, 1 - CPU%
        '''
        assert option < len(self.utime_keys)
        show_stat = show_stat if show_stat else option == BY_TIME
        self.sel_option = option
        df_set = self.get_result_df(self.sel_host)
        plot_cfg = self.PlotConfig(self.f_plot_by_eid, 'Number Parallel', self.utime_keys[option][3])
        for pair in group_eid_pairs:
            try:
                mydf = df_set if self.filters is None else self.filters(df_set)
                self.plotor.generate_stacked_group_graph(mydf, [pair[0], pair[1]], plot_cfg, pos_list, True)
            except Exception as ex:
                print("WARN, compare_group_by_eids caught exception for pair %s, ignored: %s" % (pair, ex))
                traceback.print_exc(file=sys.stdout)
        if show_stat:
            self.get_stats(None, list(itertools.chain(*group_eid_pairs)), True)
        else:
            self.show_separator()

    def plot_group_eids(self, eid_list=None, pos_list=None, option=0, show_stat=None):
        ''' eid_list = [0, 1, 2, 3] - a graph per eid, None means all
            option = 0 - wall clock, 1 - CPU%
        '''
        assert option < len(self.utime_keys)
        show_stat = show_stat if show_stat else option == BY_TIME
        self.sel_option = option
        df_set = self.get_result_df(self.sel_host)
        plot_cfg = self.PlotConfig(func=self.f_plot_by_eid, x_label='Number Parallel', y_label=self.utime_keys[option][3])
        for eid in eid_list if eid_list else df_set.eid.unique():
            try:
                self.plotor.generate_group_graph(df_set, eid, plot_cfg, pos_list, True)
            except Exception as ex:
                print("WARN, plot_group_eids caught exception for eid %d, ignored: %s" % (eid, ex))
                traceback.print_exc(file=sys.stdout)
        if show_stat:
            self.get_stats(None, eid_list, True)
        else:
            self.show_separator()

    def view_3D_by_eids(self, eid_list, pos_list=None, option=0):
        ''' NOT WORK
        my_test_name = "variant density" , NOT work
            group_eid_list = [(0, 1), (2, 3)], high value first. a graph per each (eid1, eid2)
            option = 0 - wall clock, 1 - CPU%
        '''
        assert option < len(self.utime_keys)
        self.sel_option = option
        df_set = self.get_result_df(self.sel_host)
        plot_cfg = self.PlotConfig(self.f_plot3D_by_eid, 'Number Parallel', self.utime_keys[option][3])
        try:
            plot_list = self.plotor.generate_3D_group_graph(df_set, eid_list, plot_cfg, pos_list, True)
        except Exception as ex:
            print("WARN, view_3D_by_eids caught exception for list %s: %s" % (eid_list, ex))
            traceback.print_exc(file=sys.stdout)

    def get_stats(self, gid_list=None, eid_list=None, show_stat=False):
        ''' xid = None means all'''
        if not gid_list:
            gid_list = range(1, len(VariantDensityResult.GIDLookup) + 1)
        if not eid_list:
            eid_list = range(len(VariantDensityResult.EID_Lookup))
        if self.df_stat_mean is None:
            self._build_stats_mean_df(self.sel_host)
        df = self.df_stat_mean
        # df['run name'] = df.apply(lambda r: "%s and %s mpirun" % (self.f_eid_name(r.eid), self.f_gid_name(r.gid)), axis=1)
        df = df[df['eid'].isin(eid_list) & df['gid'].isin(gid_list) & (df['num_pos'] == 1000)].sort_values(['density', 'gid'])
        stat_df = df.set_index(df.apply(lambda r: "%s and %s mpirun" % (self.f_eid_name(r.eid), self.f_gid_name(r.gid)), axis=1))
        stat_df = stat_df[['density', 'num_pos', "r/s", "rkB/s", "avgrq-sz", "avgqu-sz", "r_await", "%util"]]
        if show_stat:
            html_text = "<style>tr:nth-of-type(even) {background-color:#FAFAD2;}</style>" + tabulate.tabulate(stat_df, headers=stat_df.columns, tablefmt='html')
            display(HTML(html_text))

class VariantDensityResult(VariantDensity):
    ''' for 1000 geno '''
    filters = lambda self, df: df[(df['gid'] != 5)]           # no 88

    EID_Lookup = ['HDD with Refblock', 'SSD with Refblock', 'HDD with Varonly', 'SSD with Varonly']
    GIDLookup = ['1', '8', '16', '44', '88']
    def __init__(self, proj_name='variant-density', sel_hosts=None, root_path=None):
        super(VariantDensityResult, self).__init__(proj_name, root_path, self.GIDLookup, self.EID_Lookup)
        # 2-26 <=> 2-28 are fine; 2-27[miss 1 test]<=>2-29 is fine;
        self.sel_host = sel_hosts

class IllmnResult(VariantDensity):
    ''' for illmn '''

    EID_List = {'x': ['HDD 10000', 'SSD 10000'],    # print-calls
    'c' : ['HDD 15000', 'SSD 15000'],               # y - print-calls, c - --print-AC
    '_' : ['HDD 5000', 'SSD 5000', 'HDD 10000', 'SSD 10000', 'HDD 15000', 'SSD 15000']} # old pos pick algorithm
    GIDLookup = ['1', '4', '8', '16']           # c, y are not mpirun, mpitrun

    def __init__(self, pref, proj_name='illmn', sel_hosts=None, root_path=None):
        eid_Lookup = self.EID_List[pref] if pref in self.EID_List else self.EID_List['_']
        super(IllmnResult, self).__init__(proj_name, root_path, self.GIDLookup, eid_Lookup)
        # 2-26 <=> 2-28 are fine; 2-27[miss 1 test]<=>2-29 is fine;
        self.sel_host = sel_hosts
        self.EID_Lookup = eid_Lookup
        self.NSample2EID = [int(x.split()[1]) for x in eid_Lookup]

    def compare_by_num_sample(self, gid_list=None, num_samples=None, pos_list=None, option=0, show_stat=None):
        ''' gid_list - list of parallel [1, 8, 16]
            num_samples - [5000, 10000, 15000], eid%2=0 => HDD eid%2=1 => SSD
            option = 0 - wall clock, 1 - CPU%
        '''
        assert option < len(self.utime_keys)
        show_stat = show_stat if show_stat is not None else option == BY_TIME
        self.sel_option = option
        df_set = self.get_result_df(self.sel_host)
        plot_cfg = self.PlotConfig(self.f_plot_by_n_samples, 'Number Samples', self.utime_keys[option][3])
        try:
            eid_list = [i for i, x in enumerate(self.NSample2EID) if x in num_samples] if num_samples \
             else [n for n in range(len(self.EID_Lookup))]
            for gid in gid_list:
                rid_list = [(gid, eid) for eid in eid_list]
                # print('compare_by_num_sample: num_samples=', num_samples, ', rid_list=', rid_list, ', eid_list=', eid_list, ', gid_list=', gid_list)
                self.plotor.generate_stacked_group_graph(df_set, rid_list, plot_cfg, pos_list, True)
                if show_stat:
                    self.get_stats([gid], eid_list, True)
                else:
                    self.show_separator()
        except Exception as ex:
            print("WARN, compare_by_num_sample caught exception for gid=%s, num_samples=%s, ignored: %s" % \
            (gid_list, num_samples, ex))
            traceback.print_exc(file=sys.stdout)
