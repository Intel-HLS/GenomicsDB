#! /usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import Axes3D

# from matplotlib import (cm, pyplot as plt, mlab)

class PlotGenerator():
    # ref https://material.io/guidelines/style/color.html#color-color-palette
    names = ['Sparse', 'intermediate', 'Dense']
#    colors = [['#1A237E', '#BF360C', '#1B5E20'], ['#00796B', '#FF6F00', '#33691E']]
    colors = [['#311B92', '#FF6F00', '#33691E'], ['#00796B', '#880E4F', '#827717']] #light, dark
    colors_3d = [
    [['#311B92', '#7986CB', '#3F51B5'], ['#E65100', '#FBC02D', '#FF8F00'], ['#33691E', '#B2FF59', '#7CB342']], 
    [['#40C4FF', '#01579B', '#00B0FF'], ['#FFAB40', '#DD2C00', '#FF6E40'], ['#81C784', '#1B5E20', '#00E676']]]     #light, dark

    # top_colors = [['#304FFE', '#FF3D00', '#00C853'], ['#00BFA5', '#FFAB00', '#64DD17']]
    DefaultFontSize = 14
    fig_w, fig_h = 15, 6
    width = 0.1
    fontP = FontProperties()
    fontP.set_size('small')

    def __init__(self):
        mpl.style.use('default')
        self.xticknames = None                 # [1, 10, 100, 1000]
        self.seg_sizes = None                  # [500, 1000, 3000, 6000, 9000]
        self.alphas = None                      # [1.0, 0.4, 0.8, 0.2, 0.6] ?
        self.ind = None
        self.rects = []
        self.plot_names = []

    @staticmethod
    def calc_alpha(num_vals):
        step = 1.0/num_vals
        a = [i for i in np.arange(step, 1.0 + step, step)]
        a.reverse()
        for i in range(1, len(a) - 1, 2):
            tmp = a[i]
            a[i] = a[i+1]
            a[i+1] = tmp
        return a

    @staticmethod
    def autolabel(ax, myrects):
        """  Attach a text label above each bar displaying its height  """
        for i, rect in enumerate(myrects):
            height = rect.get_height()
            val = '%d' % int(height) if height > 2 else '%.2f' % height
            ax.text(rect.get_x() + rect.get_width()/2., 1.00*height, val, ha='center', va='bottom')

    def init_plot(self, df, xaxe):
        if self.xticknames is None:
            self.xticknames = xaxe                         # pos =[1, 10, 100, 1000] or...
            self.seg_sizes = df.seg_size.unique()            # [500, 1000, 3000, 6000, 9000]
            if len(self.seg_sizes) == 0:
                print('!!!! 0 self.seg_sizes')
                print(df)
                self.seg_sizes = [1000, 5000, 10000]
            self.alphas = self.calc_alpha(len(self.seg_sizes))  # [1.0, 0.4, 0.8, 0.2, 0.6] ?
            self.ind = np.arange(len(self.xticknames))
            self.fig_w = len(self.xticknames) * 4
        self.plot_names.clear()
        self.rects.clear()
        plt.close('all')
        _, ax = plt.subplots(1, 1, figsize=(self.fig_w, self.fig_h))
        ax.yaxis.grid(True)               #ax.grid(True) for both
        return ax

    def calc_bars(self, df, color_pallet, exp_id=''):
        # print(major_names, minor_names, my_colors, my_alphas)
        for i, density in enumerate(self.names):
            for j, seg_size in enumerate(self.seg_sizes):
                df1 = df[(df.density == i) & (df.seg_size == seg_size)]
                # print(df1.describe())
                self.plot_names.append("%s SS%s %s" % (density, seg_size, exp_id))
                ind_seg = self.ind + ((j + (i*len(self.seg_sizes))) * self.width)
                # removed std of test WC yerr=df1.wc_std
                rect = plt.bar(ind_seg, df1.wc_mean, width=self.width, align='center',\
                                color=color_pallet[i], alpha=self.alphas[j], yerr=df1.wc_std)
                self.rects.append(rect)
        # print("calc_bar, #df=%d, #rect=%d" % (len(df.index), len(self.rects)))

    def make_plot(self, ax, title, x_label, y_label):                  #, title, x_title='Number of Positions'):
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title, fontsize=self.DefaultFontSize)
        ax.set_xticks(self.ind + self.width * (len(self.seg_sizes)) * len(self.names) / 2)
        ax.set_xticklabels(self.xticknames)
        trects_0 = tuple([r[0] for r in self.rects])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        ax.legend(trects_0, tuple(self.plot_names), prop=self.fontP, loc='center left', bbox_to_anchor=(1, 0.5))
        _ = [self.autolabel(ax, r) for i, r in enumerate(self.rects) if i % 3 == 1]

    def generate_group_graph(self, df, run_id, plot_cfg, pos_list=None, show_graph=False):
        ''' run_id: either gid or eid
        '''
        assert df is not None
        my_pos_list = pos_list if pos_list else df.pos.unique()
        for i, pos in enumerate(my_pos_list):
            df_list, title, xticks = plot_cfg.func(df, pos, [run_id])
            assert len(df_list) == 1
            ax = self.init_plot(df_list[0], xticks)
            self.calc_bars(df_list[0], self.colors[i % 2])
            self.make_plot(ax, title, plot_cfg.x_label, plot_cfg.y_label)
            if show_graph:
                plt.show()

    def generate_stacked_group_graph(self, df, runid_list, plot_cfg, pos_list, show_graph=False):
        assert df is not None
        my_pos_list = pos_list if pos_list else df.pos.unique()
        for pos in my_pos_list:
            df_list, title, xticks = plot_cfg.func(df, pos, runid_list)
            assert len(df_list) == 2
            ax = self.init_plot(df_list[0], xticks)
            for i, mydf in enumerate(df_list):
                self.calc_bars(mydf, self.colors[i])
            self.make_plot(ax, title, plot_cfg.x_label, plot_cfg.y_label)
            if show_graph:
                plt.show()

    def generate_3D_group_graph(self, df, runid_list, plot_cfg, pos_list, show_graph=False):
        assert df is not None
        my_pos_list = pos_list if pos_list else df.pos.unique()
        df_list, title, xticks, zticks = plot_cfg.func(df, my_pos_list, runid_list)
        assert len(df_list) == len(runid_list) * len(my_pos_list)
        plt.close('all')
        fig = plt.figure(figsize=(len(xticks) * 6, 12))
        ax = fig.add_subplot(111, projection='3d')
        self.plot_names.clear()
        self.rects.clear()
        # ax.yaxis.grid(True)               #ax.grid(True) for both
        self.ind = np.arange(len(xticks))
        self.seg_sizes = df_list[0].seg_size.unique()            # [500, 1000, 3000, 6000, 9000]
        num_z = len(zticks)
        alphas = [1, 0.8]               #np.linspace(1, (11-num_z)/10, num=num_z)
        z_vals = np.linspace(num_z, 1, num=num_z)
        for ci, mydf in enumerate(df_list):
            for i, density in enumerate(self.names):
                for j, seg_size in enumerate(self.seg_sizes):
                    # print(" i=%d, j=%d, seg_size=%d" % (i, j, seg_size))
                    df1 = mydf[(mydf.density == i) & (mydf.seg_size == seg_size)]
                    if ci < len(self.colors_3d):
                        self.plot_names.append("%s SS%s" % (density, seg_size))
                    ind_seg = self.ind + ((j + (i*len(self.seg_sizes))) * self.width)
                    rect = ax.bar(ind_seg, df1.wc_mean, zs=z_vals[ci], zdir='y', width=self.width, align='center', color=self.colors_3d[ci % 2][i][j], alpha=0.85) 
                    self.rects.append(rect)
        ax.set_xlabel(plot_cfg.x_label)
        ax.set_zlabel(plot_cfg.y_label)
        ax.set_title(title, fontsize=self.DefaultFontSize)
        ax.set_xticks(self.ind + self.width * (len(self.seg_sizes)) * len(self.names) / 2)
        ax.set_xticklabels(xticks)                    # self.xticknames
        ax.set_yticks(z_vals)
        ax.set_yticklabels(zticks)
            
        trects_0 = tuple([r[0] for r in self.rects])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(trects_0, tuple(self.plot_names), prop=self.fontP, loc='center left', bbox_to_anchor=(1, 0.5))
        if show_graph:
            plt.show()

    # def generate_stacked_group_graph2(self, df, runid_list, plot_cfg, pos_list):
    #     assert df is not None
    #     plt_list = []
    #     my_pos_list = pos_list if pos_list else df.pos.unique()
    #     for pos in my_pos_list:
    #         df_list, title, xticks = plot_cfg.func(df, pos, runid_list)
    #         assert len(df_list) == 2
    #         ax = self.init_plot(df_list[0], xticks)
    #         for i, mydf in enumerate(df_list):
    #             self.calc_bars(mydf, self.colors[i])
    #         self.make_plot(ax, title, plot_cfg.x_label, plot_cfg.y_label)
    #         plt_list.append(plt)
    #         print("generate_stacked_group_graph2, #plot=", len(plt_list))
    #     return plt_list
# ############### for stat ???
# import math
# from matplotlib import (cm, pyplot as plt, mlab)

# def visualize(word, model):
#     """ visualize the input model for a particular word """
#     variance=np.array([np.diag(model.covars_[i]) for i in range(model.n_components)])
#     figures = []
#     for parm_idx in range(len(model.means_[0])):
#         xmin = int(min(model.means_[:,parm_idx]) - max(variance[:,parm_idx]))
#         xmax = int(max(model.means_[:,parm_idx]) + max(variance[:,parm_idx]))
#         fig, axs = plt.subplots(model.n_components, sharex=True, sharey=False)
#         colours = cm.rainbow(np.linspace(0, 1, model.n_components))
#         for i, (ax, colour) in enumerate(zip(axs, colours)):
#             x = np.linspace(xmin, xmax, 100)
#             mu = model.means_[i,parm_idx]
#             sigma = math.sqrt(np.diag(model.covars_[i])[parm_idx])
#             ax.plot(x, mlab.normpdf(x, mu, sigma), c=colour)
#             ax.set_title("{} feature {} hidden state #{}".format(word, parm_idx, i))

#             ax.grid(True)
#         figures.append(plt)
#     for p in figures:
#         p.show()

# visualize(my_testword, model)
