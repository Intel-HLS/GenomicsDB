#! /usr/bin/python3
# pylint: disable=undefined-loop-variable

import os
import os.path
import math
import logging
import pandas

logger = logging.getLogger()

class HistogramManager(object):
    # DIST_RANDOM = 'random'   # dense and sparse defined differently 
    # DIST_DENSE = 'dense'
    # DIST_SPARSE = 'sparse'

    def __init__(self, file_path):
        if os.path.isfile(file_path):
            self.df = pandas.read_csv(file_path, skiprows=1, skipfooter=1, header=None, engine='python')
            self.df.columns = ['start_pos', 'end_pos', 'byte_size']
            self.df['pos'] = range(len(self.df.index)) 
        else:
            logger.warning("file %s not found", file_path)

    def calc_bin_idx_pos(self, bin_num):
        ''' return a list of begin pos '''
        bin_size = math.ceil(self.df['byte_size'].sum() / bin_num)
        bin_start_list = []
        subtotal = 0
        for idx, row in self.df.iterrows():
            if idx == 0:
                bin_start_list.append(row)
                if bin_num == 1:
                    break
            subtotal += row['byte_size']
            if subtotal > bin_size:
                bin_start_list.append(row)
                subtotal = 0
                if len(bin_start_list) == bin_num:
                    break
        if subtotal:
            bin_start_list.append(row)
        assert len(bin_start_list) == bin_num, "not match: len(bin_start_list)=%d, bin_num=%d" %(len(bin_start_list), bin_num)
        return bin_start_list

    def calc_bin_begin_pos(self, bin_num):
        bin_start_list = self.calc_bin_idx_pos(bin_num)
        begin_list = [item['start_pos'].item() for item in bin_start_list]
        return begin_list

    # def densePos(self, mydf, num_pos):
    #     mid = math.sqrt(num_pos)
    #     v1 = int(mid + 0.5)
    #     v2 = math.ceil(mid)
    #     denselist = mydf.sort_values(['byte_size'], ascending=False).head(v1)
    #     pos_list = []
    #     for idx, row in denselist.iterrows():
    #         pos_list.extend(random.sample(range(row['start_pos'], row['end_pos']), v2))
    #     return pos_list[:num_pos]

    # def sparsePos(self, mydf, num_pos):
    #     start_pos = mydf['start_pos'].iloc[0]
    #     end_pos = mydf['end_pos'].iloc[-1]        
    #     gap = math.floor((end_pos - start_pos) / (num_pos))
    #     start = int(gap/2) + start_pos
    #     pos_list = [pos for pos in range(start, end_pos, gap)]
    #     return pos_list

    # def getPositions(self, dist_type, num_pos, from_to):
    #     eix = from_to[1] if from_to[1] else self.df.iloc[len(self.df)-1]['pos']
    #     sect_df = self.df[from_to[0] : eix]
    #     if dist_type == self.DIST_RANDOM:
    #         pos_list = random.sample(range(sect_df['start_pos'].iloc[0], sect_df['end_pos'].iloc[-1]), num_pos)
    #     elif dist_type == self.DIST_DENSE:
    #         pos_list = self.densePos(sect_df, num_pos)  
    #     elif dist_type == self.DIST_SPARSE: 
    #         pos_list = self.sparsePos(sect_df, num_pos)     
    #     else:
    #         raise ValueError
    #     return pos_list
