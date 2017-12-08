#! /usr/bin/python3

'''
wc -l merged_16p_10000.intpos => 145418043
merged_16p_10000.intpos last pos = 3036293733
    variants: 145418043 / 3036293733 4.79% 
number of base pair for human = 3.3x10**9  http://www.edinformatics.com/math_science/human_genome.htm
    variants: 145418043 / 3300000000 = 4.41%
'''
import os.path
import threading
import time
import pandas as pd
import shared.profiling_env as ppenv
import t_count_manager

NUM_LENGTH_LIST = [10, 100, 1000]
NUM_READ_ROW = 5000

intpos_fn = 'merged_16p_10000.intpos'
proj_name = 'illmn'

def process_main(fpath, len_list):
    mythread = threading.Thread(target=t_count_manager.count_manager, args=(len_list))
    mythread.start()
    fpath = os.path.join(ppenv.get_intpos_root(), proj_name, 'interesting_pos', fpath)
    num_rows = min(len_list.min)
    num_read_rows = 0
    while num_read_rows < NUM_READ_ROW:
        df = pd.read_csv(intpos_fn, skiprows=num_read_rows, nrows=num_rows, header=None, engine='python', delimiter=r"\s+")
        t_count_manager.add_df(df)
        num_read_rows += num_rows
        t_count_manager.num_read_rows = num_read_rows
    t_count_manager.done_read = True
    mythread.join()

if __name__ == '__main__':
    print('main tid=', threading.get_ident())
    process_main(intpos_fn, NUM_LENGTH_LIST)
    # test()
