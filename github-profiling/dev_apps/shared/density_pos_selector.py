#! /usr/bin/python3
import os
import os.path
import time
import logging
import pandas as pd
from profiling_env import get_cfg_json_ws, get_interest_pos_root, get_density_pos_input
from pygenomicsdblib.utils.constant_def import DENSITY_NAMES
logger = logging.getLogger()

class IntPosBuilder():
    df_list = {}
    def __init__(self, proj_name=None):
        self.interesting_pos_toproot = get_interest_pos_root(proj_name)

    def load_intpos_file(self, src_fn):
        ckpt_root = os.path.join(get_cfg_json_ws(), "checkpoints")
        if not os.path.isdir(ckpt_root):
            os.makedirs(ckpt_root)
        if not os.path.isfile(src_fn):
            src_fn_2 = os.path.join(self.interesting_pos_toproot, src_fn)
            assert os.path.isfile(src_fn_2)
            src_fn = src_fn_2
        base_name = os.path.basename(src_fn).split('.')
        base_name = base_name if len(base_name) == 1 else base_name[:-1]
        checkpoint_fn = os.path.join(ckpt_root, ''.join(base_name) +'.ckpt')
        df = None
        if os.path.isfile(checkpoint_fn):
            logger.debug("found checkpoint_fn=%s, will load from the file", checkpoint_fn)
            if os.path.getmtime(src_fn) < os.path.getmtime(checkpoint_fn):
                df = pd.read_pickle(checkpoint_fn)
                assert isinstance(df, pd.core.frame.DataFrame)
                logger.info("loaded from checkpoint_fn=%s, df_type=%s", checkpoint_fn, type(df))
            else:
                logger.info('source file %s has changed', src_fn)
                os.remove(checkpoint_fn)
        if df is None:
            logger.debug("not found %s, load from %s", checkpoint_fn, src_fn)
            df = self.load_from_file(src_fn)
            logger.debug("pickle to %s", checkpoint_fn)
            df.to_pickle(checkpoint_fn)
        return df
    @staticmethod
    def load_from_file(src_fn):
        ''' genopp_process_utils.py:get_variant_df(src_fn) load from chkpoint, not here '''
        start_time = time.time()
        logger.debug("load_from_file loading %s at %s", src_fn, time.ctime())
        df = pd.read_csv(src_fn, skiprows=1, skipfooter=1, header=None, engine='python', delimiter=r"\s+")
        elapse_time = time.time() - start_time
        logger.info("load_from_file completed at %s, duration = %d", time.ctime(), elapse_time)
        df.columns = ['position', 'num_intersect', 'num_ref_block_cell', 'num_spos_aligned']
        df['num_variant'] = df.num_intersect.sub(df.num_ref_block_cell, axis=0)
        df.drop(df[df.num_variant == 0].index, inplace=True)    #remove row with 0 variants
        logger.debug("INFO: load csv file %s in %d sec. %d rows has variant", src_fn, elapse_time, len(df.index))
        return df

    @staticmethod
    def _get_density_pos(df_all, region_pos, low=0, high=float("inf")):
        ''' 10/2, replaced by _get_random_density_pos
            df_all - dataframe of intersect file from KG
            num_pos - list of number pos per region
            low - lowest pos (inclusive), default 0,
            high - high pos (exclusive), default float("inf") (to the end)
            return list of num_region each as [([pos], [#var]),...,([pos], [#var])]
            example: pick 2, 3 position for dense, middle and sparse
            [ [([13416, 13379], [22, 8]), ([13416, 13379, 13503], [22, 8, 5])], [..mid reg.], [..sparse..] ]
        '''
        in_range = df_all[(df_all.position >= low) & (df_all.position < high)]
        mean_val = in_range.num_variant.mean()
        def mid_pos(df, nh, nl):
            return pd.concat([df[df['num_variant'] >= mean_val].tail(nh), df[df['num_variant'] < mean_val].head(nl)])
        sorted_df = in_range.sort_values(by=['num_variant', 'position'], ascending=[False, True])
        logger.debug("mean_val=%.2f, pos_range(%s, %s), #var_range=(%s, %s)  #df=%s", mean_val, low, high, \
                sorted_df['num_variant'].min(), sorted_df['num_variant'].max(), sorted_df.shape)
        all_region = list()
        all_region.append([sorted_df.head(n)[['position', 'num_variant']] for n in region_pos[0]])
        all_region.append([mid_pos(sorted_df, hn, ln)[['position', 'num_variant']] for hn, ln in [(int(n/2+n%2), int(n/2)) for n in region_pos[1]]])
        all_region.append([sorted_df.tail(n)[['position', 'num_variant']] for n in region_pos[2]])
        return all_region

    @staticmethod
    def _get_random_density_pos(df_all, region_pos, rseed, pool_perc=0.15):  # 
        ''' randomly select position in dense, midmean and sparse regions. 
            df_all - dataframe of intersect file from KG
            region_pos - list of number pos per region [[250, 1000], [250, 1000], [250, 1000]]
            pool_perc - percentile of region for the selection pool
            low - lowest pos (inclusive), default 0,
            high - high pos (exclusive), default float("inf") (to the end)
            return list of num_region each as [[dense_df250, dense_df1000], [mid_df250, mid_df1000], [sparse_df250, sparse_df1000]]
        '''
        sorted_df = df_all.sort_values(by=['num_variant', 'position'], ascending=[False, True])
        mean_val = sorted_df.num_variant.mean()
        dhigh = sorted_df.num_variant.max()
        dlow = max(dhigh - ((dhigh - mean_val) / 3), dhigh * (1 - pool_perc))
        slow = sorted_df.num_variant.min()
        shigh = min(dhigh * pool_perc, mean_val / 3)
        mhigh = min(mean_val + (dhigh * pool_perc), mean_val + ((dhigh - mean_val) / 3))
        mlow = max(mean_val - (dhigh * pool_perc), mean_val * 2 / 3)
        logger.debug("mean_val=%.2f, var_high=(%s, %s), var_low=(%s, %s), var_mid=(%s, %s), pos_range(%s, %s), #df=%s", \
            mean_val, dhigh, dlow, shigh, slow, mhigh, mlow, sorted_df.position.min(), sorted_df.position.max(), sorted_df.shape)
        all_region = list()
        all_region.append([sorted_df[(sorted_df.num_variant >= dlow)].sample(n, random_state=rseed)[['position', 'num_variant']] for n in region_pos[0]])
        all_region.append([sorted_df[(sorted_df.num_variant <= mhigh) & (sorted_df.num_variant >= mlow)].sample(n, random_state=rseed)[['position', 'num_variant']] for n in region_pos[1]])
        all_region.append([sorted_df[(sorted_df.num_variant <= shigh)].sample(n, random_state=rseed)[['position', 'num_variant']] for n in region_pos[2]])
        _ = [[logger.info("%s, %s", DENSITY_NAMES[i], pl.describe()) for pl in dpl] for i, dpl in enumerate(all_region)]
        return all_region
 
    def get_query_pos_list(self, mydf, region_pos, array_start_pos, num_parallel):
        ''' get pos in arrays for query cfg for dense, mid and sparse respestively for all num_pos to pick
        columns_list: [ [([[pos_4 arr1],.. ], [max_n_val, min_n_val]), ..(next,num_pick)], #dense
                [([[pos_4 arr1],.. ], [max_n_val, min_n_val]), ..(next,num_pick)],  #mid
                [([[pos_4 arr1],.. ], [max_n_val, min_n_val]), ..(next,num_pick)] ] #sparse
        '''
        def merge_sel_pos(merge_list, df):
            merge_list[0].append(df['position'].tolist())
            merge_list[1].append(df['num_variant'].min())
            merge_list[2].append(df['num_variant'].max())

        #print(region_pos)
        column_sel_list = None             #[[]] * len(region_pos) =>[[],[],[]]
        save2list = list()
        for i in range(min(num_parallel, len(array_start_pos))):
            spos = array_start_pos[i]
            epos = array_start_pos[i+1] if i+1 < len(array_start_pos) else float("inf")
            # old algorithm _get_density_pos
            sel_pos_lists = self._get_random_density_pos(mydf[(mydf.position >= spos) & (mydf.position < epos)], region_pos, num_parallel)
            save2list.extend(sel_pos_lists)

            if not column_sel_list:
                column_sel_list = [[([x['position'].tolist()], [x['num_variant'].min()], [x['num_variant'].max()]) for x in dfl] for i, dfl in enumerate(sel_pos_lists)]
            else:
                _ = [[merge_sel_pos(column_sel_list[i][j], x) for j, x in enumerate(dfl)] for i, dfl in enumerate(sel_pos_lists)]
        # print(column_sel_list)
        return column_sel_list, save2list

    def get_intpos_df(self, num_sample=None, num_parallel=None):
        fn = get_density_pos_input()
        if fn not in self.df_list:
            self.df_list[fn] = self.load_intpos_file(fn)
        return self.df_list[fn], fn

class IllmnIntPosBuilder(IntPosBuilder):
    def __init__(self):
        super(IllmnIntPosBuilder, self).__init__('illmn')
        self.df_intpos = self.__build_intpos_info(self.interesting_pos_toproot)

    @staticmethod
    def __build_intpos_info(topdir):
        sample_dirs = [(sd, sd.split('sample')[1]) for sd in os.listdir(topdir) if os.path.isdir(os.path.join(topdir, sd)) and sd.startswith('sample')]
        tlist = []
        for item in sample_dirs:     # (dirname, num_sample)
            n_sample = int(item[1][:-3])
            iposfiles = os.listdir(os.path.join(topdir, item[0]))
            iposfnlist = sorted(iposfiles, key=lambda x: int(x.split('_')[-1]))
            tlist.append([n_sample, item[0], iposfnlist])
        df_intpos = pd.DataFrame(tlist, columns=["n_sample", "dname", 'files'])
        df_intpos['combined'] = df_intpos.apply(lambda r: 'merged_%sp_%s.intpos' % (len(r['files']), r['n_sample']), axis=1)
        # logger.debug(df_intpos)
        return df_intpos

    def _build_intpos(self, num_sample):
        selected = self.df_intpos[(self.df_intpos.n_sample <= num_sample)].sort_values(by='n_sample', ascending=False).head(1)
        # assert len(selected.iloc[0]['files']) >= num_parallel
        ret_fn = os.path.join(self.interesting_pos_toproot, selected.iloc[0]['combined'])
        if os.path.isfile(ret_fn):
            print("INFO found %s for num_sample=%d" % (ret_fn, num_sample))
            return ret_fn
        try:
            with open(ret_fn, "w") as outfile:
                num_line = 0
                for fn in selected.iloc[0]['files']:
                    ipfn = os.path.join(self.interesting_pos_toproot, selected.iloc[0]['dname'], fn)
                    with open(ipfn, "r") as infile:
                        for aline in infile:
                            try:
                                _ = [int(token) for token in aline.split()]
                                num_line += 1
                                outfile.write(aline)
                            except ValueError:
                                print("INFO skip ", aline)
            print("INFO made %s contains %d lines for num_sample=%d" % (ret_fn, num_line, num_sample))
        except Exception:
            if os.path.isfile(ret_fn):
                os.remove(ret_fn)
            return None
        return ret_fn

    def get_intpos_df(self, num_sample=1000):
        ''' 
        return df for the num_sample and num_parallel
        '''
        assert num_sample
        fn = self._build_intpos(num_sample)
        if fn not in self.df_list:
            self.df_list[fn] = self.load_intpos_file(fn)
        return self.df_list[fn], fn

def get_intpos_builder(proj_name):
    return IllmnIntPosBuilder() if proj_name == 'illmn' else IntPosBuilder()
