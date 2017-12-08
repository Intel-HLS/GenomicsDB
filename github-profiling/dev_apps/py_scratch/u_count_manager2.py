#! /usr/bin/python3

'''
26100003
35321003
145418043
wc -l merged_16p_10000.intpos => 145418043
merged_16p_10000.intpos last pos = 3036293733
    variants: 145418043 / 3036293733 4.79% 
number of base pair for human = 3.3x10**9  http://www.edinformatics.com/math_science/human_genome.htm
    variants: 145418043 / 3300000000 = 4.41%

'''
import os
import os.path
import threading
import time
import platform
import pickle
import pandas as pd
import shared.profiling_env as ppenv

LOGGER_ROOT = '/path/to/logs'
is_test = "CM_TEST" in os.environ
if is_test:
    CHECK_FREQ = 1                   # a sec
    NUM_LENGTH_LIST = [10, 35]
    intpos_fn = '1000g_interesting_positions.200'
    proj_name = 'query_density'
else:
    CHECK_FREQ = 60 * 15         # 15 minutes
    NUM_LENGTH_LIST = [10000, 100000, 1000000]                      #[10, 100, 1000]
    intpos_fn = 'merged_16p_10000.intpos'
    proj_name = 'illmn'

Print2ScreenNum = 3

def counter_all(mydf):
    ''' no thread pool '''
    log_msg('++counter: mydf={}'.format(mydf.shape))
    drstl = []                   # [pos, sum_rl1, sum_rl2, sum_rl3]
    sidx = [0] * len(NUM_LENGTH_LIST)
    def show_log(ridx):         
        ''' 145418043 / 500000 = 291 '''
        return ridx > 1 if is_test else (ridx + 1) % 500000 == Print2ScreenNum
    def get_spos():
        return ";".join(["{}-{}:{}".format(NUM_LENGTH_LIST[i], drstl[ii][0], drstl[ii][i+1]) for i, ii in enumerate(sidx)])
    for ridx, row in mydf.iterrows():
        # print(ridx, row)
        cpos = row[0]
        nvar = row[3]
        if show_log(ridx):
            last_3 = drstl[Print2ScreenNum * -1:]
            log_msg('B: #drstl={}, ridx={}, sidx={}, spos={}, cpos=({},{}), last_3={}'.format(len(drstl), ridx, sidx, get_spos(), cpos, nvar, drstl[Print2ScreenNum * -1 : ]))
        stime = time.time()
        idx = min(sidx)
        while idx < len(drstl):
            for i, rlen in enumerate(NUM_LENGTH_LIST):
                if idx >= sidx[i]:
                    if drstl[idx][0] + rlen > cpos:
                        drstl[idx][1 + i] += nvar
                    else:
                        sidx[i] = idx
            idx += 1
        drstl.append([cpos] + [nvar] * len(NUM_LENGTH_LIST))              
        if show_log(ridx):
            log_msg('A: {}sec: #drstl={}, sidx={}, spos={}, last_3={}'.format((time.time() - stime), len(drstl), sidx, get_spos(), drstl[Print2ScreenNum * -1 : ]))
    return drstl

def process_all(ip_fn):
    start_time = time.time()
    fpath = os.path.join(ppenv.get_intpos_root(), proj_name, 'interesting_pos', ip_fn)
    df_all = pd.read_csv(fpath, header=None, engine='python', delimiter=r"\s+")
    log_msg('Done read file in {} sec, #main={}'.format((time.time() - start_time), df_all.shape))

    results = counter_all(df_all)
    save_fn = os.path.join(LOGGER_ROOT, "{}-u_{}.chkpnt".format(ip_fn, time.strftime('%m%d%H%M', time.localtime())))
    with open(save_fn, 'wb') as ofd:
        pickle.dump(results, ofd, protocol=pickle.HIGHEST_PROTOCOL)
    log_msg('process_all: Done in {} sec, #in_rec={}, #out_rec={}, saves to {}'.format(time.time() - start_time, df_all.shape, len(results), save_fn))

log_to = {}
def log_msg(msgstr):
    tid = threading.get_ident()
    try:
        if is_test:
            print(msgstr)
        with open(log_to[tid], 'a') as lfd:
            lfd.write('{}@{}: {} \n'.format(tid, time.strftime('%y/%m/%d %H:%M:%S', time.localtime()), msgstr))
    except Exception as ex:
        if tid not in log_to or not os.path.exists(log_to[tid]):
            init_logger()
            with open(log_to[tid], 'a') as lfd:
                lfd.write('{}@{}: {} \n'.format(tid, time.strftime('%y/%m/%d %H:%M:%S', time.localtime()), msgstr))
        else:
            print("log_msg: ignored exception {}".format(ex))             

def init_logger(isMain=False):
    tid = threading.get_ident()
    if tid in log_to and os.path.exists(log_to[tid]):
        print("Info: for tid {}, log file {} already created".format(tid, log_to[tid]))
    else:
        log_fn = os.path.join(LOGGER_ROOT, "rnv_{}-{}-{}_2.log".format('main' if isMain else 'other', proj_name, platform.node().split('.')[0][-2:])) 
        log_to[tid] = log_fn
        with open(log_to[tid], 'w') as lfd:
            lfd.write('{}@{}: Starts.. \n'.format(tid, time.strftime('%y/%m/%d %H:%M:%S', time.localtime())))
        print("Created log file {} for tid {}".format(log_to[tid], tid))

if __name__ == '__main__':
    my_tid = threading.get_ident()
    init_logger(True)
    process_all(intpos_fn)
    print('Done: logs: {}'.format(';'.join(['{} to fn {}'.format("main" if k == my_tid else "other", os.path.basename(v)) for k, v in log_to.items()])))
