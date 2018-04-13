#! /usr/bin/python3
# pylint: disable=undefined-loop-variable

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
import math
import logging
import pandas

logger = logging.getLogger()

class HistogramManager(object):

    def __init__(self, file_path):
        if os.path.isfile(file_path):
            df = pandas.read_csv(file_path, skiprows=1, skipfooter=1, header=None, engine='python')
            df.columns = ['start_pos', 'end_pos', 'byte_size']
            df['pos'] = range(len(df.index))
            self.df = df
        else:
            logger.warning("file %s not found", file_path)

    def calc_bin_idx_pos(self, bin_num, begin_pos=1, end_pos=-1):
        ''' return a list of begin pos '''
        df = self.df[self.df['end_pos'] >= begin_pos] if end_pos == -1 else self.df[(self.df['start_pos'] >= begin_pos) & (self.df['end_pos'] <= end_pos)]
        bin_size = math.ceil(df['byte_size'].sum() / bin_num)
        bin_start_list = []
        subtotal = 0
        for idx, row in df.iterrows():
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

    def calc_bin_begin_pos(self, bin_num, begin_pos=1, end_pos=-1):
        bin_start_list = self.calc_bin_idx_pos(bin_num, begin_pos, end_pos)
        begin_list = [item['start_pos'].item() for item in bin_start_list]
        return begin_list
