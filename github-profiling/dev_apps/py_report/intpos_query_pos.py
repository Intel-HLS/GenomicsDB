import sys
import os.path
import re
from collections import namedtuple
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pygenomicsdblib.utils.constant_def import DENSITY_NAMES
from shared.density_pos_selector import get_intpos_builder

# FORMAT = '%(asctime)-15s %(clientip)s %(user)-8s %(message)s'
# logging.basicConfig(format=FORMAT)
# logger = logging.getLogger()
# logger.setLevel(logging.DEBUG)

DefaultConfigDir = "C:\\Users\\mrutarx\\myprojects\\GenomicsDBPerfTest\\run_tests"
PathPatten = namedtuple('PathPatten', ["all_dirs", "ignore_tdb_dirs", "cfg_fname"])
#ProjectPref = 'i'

pen_sizes = [0.1, 8, 40]                        # fine, middle, thick
# color for interesting positions and selected positions
color_ip = '#EFEBE9'                            # interesting positions
color_mean = '#455A64'
#g, r, b
color_nv = [['#76FF03', '#FF3D00', '#2196F3'], ['#D4E157', '#FFAB40', '#18FFFF']]    
    # color for intpos-selpos color = ['#311B92', '#FF6F00', '#33691E'] g,r,b <=> dense, mid, sparse
# color for selected positions
color = ['#2E7D32', '#F44336', '#304FFE']       # color for selpos = ['#311B92', '#FF6F00', '#33691E'] g,r,b <=> dense, mid, sparse
color_edge ='#777777'
color_line = "#cccccc"
MeanLW = 5
markers = ['^', 'v', 'o', '+']

class ViewPositions():
    IllmnTestsPatterns = {
        'i': PathPatten(all_dirs="^i-[12][0-5]_16-[14]_[sh][sd]d_.*_1{0,1}[50]000$", ignore_tdb_dirs="^i-[12][0-5]_16-[14]_ssd_.*_1{0,1}[50]000$", cfg_fname="^[14]-qc_.*-1{0,1}[50]000.json$"), 
        'a': PathPatten(all_dirs="^a-[12][0-7]_16-[14]_[sh][sd]d_.*_1{0,1}[0-9][05]00$", ignore_tdb_dirs="^a-[12][0-7]_16-[14]_hdd_.*_1{0,1}[0-9][05]00$", cfg_fname="^[14]-qc_.*-1{0,1}[0-9][50]00.json$"),
        'x': PathPatten(all_dirs="^x-[1-4][01]_16-[148]6{0,1}_[sh][sd]d_.*_10000$", ignore_tdb_dirs="^x-[1-4][01]_16-[148]6{0,1}_ssd_.*_10000$", cfg_fname="^[148]6{0,1}-qc_.*-10000.json$")}
    density_lu = lambda x: 0 if x == 'dense' else 1 if x == 'meanmid' else 2 if x == 'sparse' else - 1

    def __init__(self, proj_name="illmn", jsoncfg_root=None):
        self.cfg_root_dir = jsoncfg_root if jsoncfg_root else DefaultConfigDir
        assert os.path.isdir(self.cfg_root_dir), '%s not found' % self.cfg_root_dir
        self.project_name = proj_name
        self.intpos_root = os.path.join(DefaultIntPosRoot, self.project_name, 'interesting_pos')

        self.my_density_names = [d[0].upper() for d in DENSITY_NAMES]
        self.my_density_names.reverse()
        # TODO: temporary. hardcoded
        self.intpos_subdir = 'sample10000_28'
        self.cfg_file_list = None
        self.target_cfg = None
        self.intpos_dict = {}
        pd.set_option('display.width', 160)

    @staticmethod
    def __find_cfg_files(myroot_dir, qinclude_dirs, qinclude_files):
        qcfg_files = []
        for root, dirs, files in os.walk(myroot_dir, topdown=True):
            dirs[:] = [d for d in dirs if re.match(qinclude_dirs, d)]
            # print(files)
            files = [(os.path.basename(root), f) for f in files if re.match(qinclude_files, f)]
            qcfg_files.extend(files)
        # print("#=", len(cfg_files))
        return qcfg_files

    def build_qjsoncfg_df(self, pref, proj_name=None):
        assert pref and pref in self.IllmnTestsPatterns

        if proj_name:
            self.project_name = proj_name
        incpatn = self.IllmnTestsPatterns[pref]
        rpath = os.path.join(self.cfg_root_dir, self.project_name)
        assert os.path.isdir(rpath), rpath
        query_info = self.__find_cfg_files(rpath, incpatn.ignore_tdb_dirs, incpatn.cfg_fname)
        df_cfg = pd.DataFrame([[int(di[0][2]), int(di[0][3]), int(di[1].split('-')[-1]), int(di[-2])+1, int(di[-1]),\
        ViewPositions.density_lu(fi[0]), int(fi[1]), dn, fn] \
        for di, fi, dn, fn in [(dn.split('_'), fn.split('_')[1:3], dn, fn) for dn, fn in query_info]], \
        columns=['gid', 'eid', 'parallel', 'num_frag', 'frag_size', 'density', 'num_pos', 'dname', 'fname'])
        def get_query_pos(row):
            fn = os.path.join(self.cfg_root_dir, proj_name, row['dname'], row['fname'])
            with open(fn, 'r') as ifd:
                qc = json.load(ifd)
            return qc['query_column_ranges']     # dataframe canot have n-dim np.array
        df_cfg['pos_list'] = df_cfg.apply(get_query_pos, axis=1)
        return df_cfg

    def get_pos_num_var(self, dfpos, num_pos, part_num):
        '''g[1].pos_list size=3; each has N list, N = parallel
        pos_lst = [list(itertools.chain(*pp)) for pp in g[1].pos_list] 
        return: num_parallel : [dense_df, mid_df, sparse_df] df from df_intpos[num_part] '''
        if part_num not in self.intpos_dict:     #get intpos for the partition if not have
            self.get_intpos_list([part_num])
        df_intpos = self.intpos_dict[part_num]
        return {g[0] : [df_intpos[(df_intpos.position.isin(dpl[part_num]))].drop_duplicates() for dpl in g[1].pos_list] for g in dfpos[(dfpos.num_pos == num_pos)].groupby(['parallel'])}

    def get_intpos_list(self, np_list):
        ''' np_list likes [0, 3] '''
        for p_num in np_list:
            if p_num not in self.intpos_dict:
                fn = os.path.join(self.intpos_root, self.intpos_subdir, "interesting_pos_%d" % p_num)
                assert os.path.isfile(fn), '%s not found' % fn
                self.intpos_dict[p_num] = get_intpos_builder('illmn').load_from_file(fn)

    def plot_on_intpos(self, dfjsonpos, num_pos, part_num, plot_np_list=[1]):
        ''' dfjsonpos - df of pos from query json 
            num_part - partition number, such as 0, 3, or ..
            num_pos - 250 or 1000
            TODO: 3) animation?
            num_part = 1 or 2... 
        '''
        lnpl = len(plot_np_list)
        assert plot_np_list and lnpl > 0
        myplot_np_list = plot_np_list if lnpl == 1 else plot_np_list[:2]
        if lnpl > 2:
            print('only show parallel %s' % myplot_np_list)

        dfnp_dict = self.get_pos_num_var(dfjsonpos, num_pos, part_num) # for all parellel, 1, 4, ..

        mydf_intpos = self.intpos_dict[part_num]
        plt.close()
        _, myax = plt.subplots(figsize=(20, 8))
        myax.scatter(x=mydf_intpos.position, y=mydf_intpos.num_variant, s=pen_sizes[0], c=color_ip, label='interesting position', zorder=1)
        myax.axhline(y=mydf_intpos.num_variant.mean(), c=color_mean, linewidth=MeanLW, ls='-.', label='mean value', zorder=2)

        self.myax = myax
        self.dfnp_dict = dfnp_dict
        for np in myplot_np_list:
            mycolors = color_nv[np % 2]
            for i, dfnv in enumerate(self.dfnp_dict[np]):  #densities
                self.myax.scatter(x=dfnv.position, y=dfnv.num_variant, marker=markers[i], s=pen_sizes[2], c=mycolors[i], edgecolor=color_edge, label=DENSITY_NAMES[i], zorder=10)

        #fig.tight_layout()
        myax.set_xlabel('Positions')
        myax.set_ylabel('Number of Variants')
        myax.set_title("Distribution of %d Seleced Positions in Partition %d" % (num_pos, part_num))

        myax.legend(bbox_to_anchor=('lower center'), loc=3, ncol=5)
        plt.show()

    def set_y_ticks(self, myax):
        myax.set_yticks(range(len(self.my_density_names)))
        myax.set_yticklabels(self.my_density_names)

    def plot_sel_pos_run(self, dfpos, num_pos, partition_board):
        ''' plot a graph per run '''
        num_axis = len(dfpos.parallel.unique())
        _, ax_list = plt.subplots(nrows=num_axis, ncols=1, sharex=True, sharey=True, figsize=(20, num_axis * 1))  #
        for i, g in enumerate(dfpos[(dfpos.num_pos == num_pos)].groupby(['parallel'])):
            myaxall = ax_list[i]
            # draw vertical partition lines
            for xc in partition_board:
                myaxall.axvline(x=xc, c=color_line, linewidth=1, ls='--')
            for j, dpp in enumerate(g[1]['pos_list']):
                mycolor = color[j]
                for x in dpp:       #density per partition
                    y = np.zeros_like(x) + len(self.my_density_names) - j - 1   # dense is 0, reverse it
                    myaxall.scatter(x, y, c=mycolor, s=pen_sizes[0])
            if i == 0:
                myaxall.set_title("Selected Positions for Parallel 1, 4, 8 and 16")
                self.set_y_ticks(myaxall)
                myaxall.set_xticks(partition_board)
                myaxall.set_xticklabels(range(2, 17))
        ax_list[-1].set_xlabel('Positions and Partitions')
        plt.show()

    def plot_sel_pos_partition(self, dfpos, num_pos):
        ''' plot a graph per partition '''
        for g in dfpos[(dfpos.num_pos == num_pos)].groupby(['parallel']):
            hsize = g[0] * 1 if g[0] == 1 else g[0] * 1.5
            fig, ax_pp_list = plt.subplots(g[0], sharey=True, figsize=(20, hsize))  # sharex=True, sharey=True
            fig.subplots_adjust(hspace=.6)
            for j, dpp in enumerate(g[1]['pos_list']):
                mycolor = color[j]
                for k, x in enumerate(dpp):     # for each partition
                    y = np.zeros_like(x) + len(self.my_density_names) - j - 1
                    myaxpp = ax_pp_list if g[0] == 1 else ax_pp_list[k]
                    if k == 0:
                        myaxpp.set_title("%d Partition" % g[0])
                        self.set_y_ticks(myaxpp)
        #             print("  partition {}, den={}, {} - {}".format(g[0], id, min(x), max(x)))
                    myaxpp.scatter(x, y, c=mycolor, s=pen_sizes[1], marker=markers[3])
            lastax = ax_pp_list if g[0] == 1 else ax_pp_list[-1]
            lastax.set_xlabel("Positions in a Partition")
        plt.show()
        
    # def plot_positions(target_cfg, pref, intpos_handler):
    # pos_viewer = ViewPositions()
    # df_qjson = pos_viewer.build_qjsoncfg_df(pref, target_cfg.test_name)
    # array_start_pos = target_cfg.array_start_pos

    # num_parallel = df_qjson.parallel.max()
    # intpos_dict = {}
    # for frag in target_cfg.num_sample_list:
    #     for i in range(frag[0] + 1):
    #         dbsize = (i + 1) * frag[1]
    #         df_intpos, fn = intpos_handler.get_intpos_df(dbsize, num_parallel)
    #         if fn not in intpos_dict:
    #             intpos_dict[fn] = (df_intpos, dbsize)
    # print(intpos_dict)
    
#    make df, add ind in df_cfg
# if __name__ == '__main__':
#     plot_positions(MyIllmnConfig(ProjectPref), 'i', get_intpos_builder('illmn'))
