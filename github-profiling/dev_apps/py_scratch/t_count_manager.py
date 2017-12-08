#! /usr/bin/python3

'''
wc -l merged_16p_10000.intpos => 145418043
merged_16p_10000.intpos last pos = 3036293733
    variants: 145418043 / 3036293733 4.79% 
number of base pair for human = 3.3x10**9  http://www.edinformatics.com/math_science/human_genome.htm
    variants: 145418043 / 3300000000 = 4.41%

run program: about 10000 per sec, so it will take 5 days.

'''
import sys
import os
import os.path
import threading
import time
import platform
import concurrent.futures
import psutil
import pandas as pd
import shared.profiling_env as ppenv

LOGGER_ROOT = '/path/to/logs'
str_test = "CM_TEST"                       # export CM_TEST = 1, 28 has it
if str_test in os.environ:
    CHECK_FREQ = 1                   # a sec
    NUM_LENGTH_LIST = [10, 35]
    intpos_fn = '1000g_interesting_positions.200'
    proj_name = 'query_density'
else:
    CHECK_FREQ = 60 * 15         # 15 minutes
    NUM_LENGTH_LIST = [10000, 100000, 1000000]                      #[10, 100, 1000]
    intpos_fn = 'merged_16p_10000.intpos'
    proj_name = 'illmn'

# NUM_READ_ROW = float("inf")
bcp = len(sys.argv) > 1
MaxThread = psutil.cpu_count() * 2

my_text = ""
def set_text(txt):
    global my_text
    my_text = txt

def get_text():
    return my_text, 503

def counting(spos, rlidx):
    ''' 1st = pos; 3 fields per an interval: last pos, sum of the num_var within a window, # of rec '''
    pos_val = df_main.loc[spos, 0]
    to_pos = pos_val + NUM_LENGTH_LIST[rlidx]
    with alock:
        df = df_main[(df_main[0] >= pos_val) & (df_main[0] < pos_val + NUM_LENGTH_LIST[rlidx])]
    if spos not in df_rstl.index:
        df_rstl.loc[spos, 0] = df.loc[spos, 0]
    ncol = rlidx * 3 + 1
    df_rstl.loc[spos, ncol : ncol + 2] = [to_pos, df[3].sum(), df.shape[0]]
    # log_msg('!!! counting: ridx={}, #df_rslt={} pos={} - {}, sum={}, #={}'.format(spos, len(df_rstl.index), pos_val, to_pos, sum_val, cnt_val))

def counting2(ridx, rlidx):
    ''' 1st = pos; 3 fields per an interval: last pos, sum of the num_var within a window, # of rec '''
    intval = NUM_LENGTH_LIST[rlidx]
    mydf = df_main
    pos_val = mydf.loc[ridx, 0]
    pos_val_to = pos_val + intval
    sum_val = 0
    cnt_val = 0
    
    # log_msg('+1+ counting2: ridx={}, rlidx={}, pos={} - {}'.format(ridx, rlidx, pos_val, pos_val_to))
    for i in range(intval):
        ri = ridx + i
        if mydf.loc[ri, 0] < pos_val_to:
            sum_val += mydf.loc[ri, 3]
            cnt_val += 1
        else:
            break 

    if ridx not in df_rstl.index:
        df_rstl.loc[ridx, 0] = mydf.loc[ridx, 0]
    ncol = rlidx * 3 + 1
    df_rstl.loc[ridx, ncol : ncol + 2] = [pos_val_to, sum_val, cnt_val]
    if ridx % 1000 == 1:
        log_msg('+2+ counting2: ridx={}, #df_rslt={} pos={} - {}, sum={}, #={}'.format(ridx, len(df_rstl.index), pos_val, pos_val_to, sum_val, cnt_val))

df_main = pd.DataFrame()
alock = threading.RLock()
done_read = False
df_rstl = pd.DataFrame(columns=list(range(0, len(NUM_LENGTH_LIST) * 3 + 1)), dtype=float)
    
def count_manager():
    ''' no thread pool '''
    init_logger()
    count_ctl = {x: [0, []] for x in NUM_LENGTH_LIST}
    time.sleep(0.001)
    while df_main.empty:
        time.sleep(0.0001)
    last_n_read = 0
    with alock:
        num_rows = df_main.shape[0]
    log_msg('! count_manager: df_main is ready #=%d' % num_rows)
    sys.stdout.flush()
    while not done_read or last_n_read < num_rows:
        if last_n_read < num_rows:
            for i, rlen in enumerate(NUM_LENGTH_LIST):
                last_row_idx = count_ctl[rlen][0]
                if last_row_idx <= (num_rows - rlen):
                    count_ctl[rlen][0] = num_rows - rlen + 1
                    # counting2(last_row_idx, i)
                    for ridx in range(last_row_idx, count_ctl[rlen][0]):
                        _ = counting(ridx, i) if bcp else counting2(ridx, i)
                    log_msg('count_manager: num_rows=%d, rlen=%d, last_row_idx=%d, #item=%d' % (num_rows, rlen, last_row_idx, count_ctl[rlen][0]))
                    sys.stdout.flush()
            last_n_read = num_rows
        with alock:
            num_rows = len(df_main.index)
    log_msg('Done, count_manager, main shape={}, result shape={}'.format(df_main.shape, df_rstl.shape))
    sys.stdout.flush()

def count_manager_tp():
    ''' with thread pool, slower than not use thread pool '''
    count_ctl = {x: [0, []] for x in NUM_LENGTH_LIST}
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=MaxThread)
    def check_done():
        ret = False
        for rlen in NUM_LENGTH_LIST:
            q = count_ctl[rlen][1]
            if q:
                dl = [i for i, fs in enumerate(q) if fs.done()]
                if dl:
                    ex_dict = {}
                    for qi in dl:
                        ex = q[qi].exception()
                        if ex:
                            ex_dict[ex] = ex_dict[ex] + 1 if ex in ex_dict else 1
                    ex_total = 0
                    for kex, num in ex_dict.items():
                        log_msg(' $$$ #exception {} is {}'.format(kex, num))
                        ex_total += num                      
                    log_msg(' $$$ check_done: for %d, was %d, %d are done, %d exceptions' % (rlen, len(q), len(dl), ex_total))
                    sys.stdout.flush()
                    dl.reverse()
                    _ = [q.pop(i) for i in dl]
            ret = ret or not q
        return ret
    time.sleep(0.001)
    while df_main.empty:
        time.sleep(0.001)
    last_n_read = 0
    with alock:
        num_rows = len(df_main.index)
    log_msg('! count_manager: df_main is ready #=%d' % num_rows)
    sys.stdout.flush()
    while not done_read or last_n_read < num_rows:
        if last_n_read < num_rows:
            for i, rlen in enumerate(NUM_LENGTH_LIST):
                last_row_idx = count_ctl[rlen][0]
                if last_row_idx <= (num_rows - rlen):
                    count_ctl[rlen][0] = num_rows - rlen + 1
                    # counting2(last_row_idx, i)
                    for ridx in range(last_row_idx, count_ctl[rlen][0]):
                        # log_msg('count_manager: num_rows=%d, rlen=%d, last_row_idx=%d, #item=%d, ridx=%d, #q=%d' % (num_rows, rlen, last_row_idx, count_ctl[rlen][0], ridx, len(count_ctl[rlen][1])))                        
                        if bcp:
                            fs = executor.submit(counting, ridx, i)
                        else:
                            fs = executor.submit(counting2, ridx, i)
                        count_ctl[rlen][1].append(fs)
                    log_msg('count_manager: num_rows=%d, rlen=%d, last_row_idx=%d, #item=%d, #q=%d' % (num_rows, rlen, last_row_idx, count_ctl[rlen][0], len(count_ctl[rlen][1])))
                    sys.stdout.flush()
            last_n_read = num_rows
        else:
            check_done()
        with alock:
            num_rows = len(df_main.index)
    log_msg('count_manager: finish submitting, df_main shape={}, df_rstl shape={}'.format(df_main.shape, df_rstl.shape))
    if not check_done():
        time.sleep(30)
        check_done()
    for x in NUM_LENGTH_LIST:
        concurrent.futures.wait(count_ctl[x][1])
    log_msg('Done, count_manager, df_main shape={}, df_rstl shape={}'.format(df_main.shape, df_rstl.shape))
    sys.stdout.flush()

def read_file(fpath, num_rows):
    global done_read, df_main
    while True:
        df = pd.read_csv(fpath, skiprows=len(df_main.index), nrows=num_rows, header=None, engine='python', delimiter=r"\s+")
        df_temp = df_main.append(df, ignore_index=True)
        log_msg('process_main: #read={}, #main={}, #temp={}'.format(df.shape, df_main.shape, df_temp.shape))
        with alock:
            df_main = df_temp
    done_read = True
    log_msg('Done read file, #add={}, #main={}'.format(df.shape, df_main.shape))

def read_chunk(fpath, num_rows):
    global done_read, df_main
    chunker = pd.read_csv(fpath, chunksize=num_rows, header=None, delimiter=r"\s+")
    start_time = time.time()
    for i, rows in enumerate(chunker):
        log_msg("i={}, rows={}".format(i, rows.shape))
        df_temp = df_main.append(rows, ignore_index=True)
        with alock:
            df_main = df_temp
        if i % 10 == 0:
            log_msg(' {}> {} sec, add {}, #main={}'.format(i, (time.time() - start_time), rows.shape, df_main.shape))
            sys.stdout.flush()
        start_time = time.time()
        time.sleep(0.0001)
    chunker.close()
    done_read = True
    log_msg('Done read_chunk #rec={}'.format(df_main.shape))

# def read_chunk_2(fpath, num_rows):
#     global done_read
#     chunker = pd.read_csv(fpath, chunksize=num_rows, header=None, delimiter=r"\s+")
#     def add_chunk(chk):
#         global bDone, df_main
#         df_temp = df_main.append(chk, ignore_index=True)
#         with alock:
#             df_main = df_temp
#         bDone = len(df_main.index) >= NUM_READ_ROW
#         log_msg('add_chunk: #rows=%d, #chk=%d, type=%s, bDone=%s' % (len(df_main.index), len(chk.index), type(chk), bDone))
#     [add_chunk(ck) for ck in chunker if not bDone]
#     log_msg('read_chunk: DONE read, #rows=%d, #df_main=%d' % (len(df_main.index), len(df_main.index)))
#     done_read = True

def log_msg_selected_pos(num_pos=10):
    for rl in range(1, len(df_rstl.columns), 3):
        df_sorted = df_rstl.sort_values([rl+1], axis=0, inplace=False, ascending=False)
        log_msg('-- Dense rl={}, range_len={}, #df={} --'.format(rl, NUM_LENGTH_LIST[int(rl/3)], df_sorted.shape[0]))
        log_msg(df_sorted.head(num_pos)[[0, rl, rl+1, rl+2]])
        sys.stdout.flush()
        df_sorted = df_sorted[pd.notnull(df_sorted[rl+1])]
        log_msg('-- Sparse --')
        log_msg(df_sorted.tail(num_pos)[[0, rl, rl+1, rl+2]])
        sys.stdout.flush()

def process_main(ip_fn, len_list):
    start_time = time.time()
    mythread = threading.Thread(target=count_manager)       # , args=(len_list)
    mythread.start()
    fpath = os.path.join(ppenv.get_intpos_root(), proj_name, 'interesting_pos', ip_fn)
    num_rows = max(len_list)
    # read_file(fpath, num_rows)
    read_chunk(fpath, num_rows)
    while mythread.is_alive():
        mythread.join(CHECK_FREQ)
        log_msg("check: df_rstl {}".format(df_rstl.shape))
    save_fn = os.path.join(ppenv.get_intpos_root(), proj_name, 'interesting_pos', "%s.saved" % ip_fn)
    df_rstl.to_pickle(save_fn)
    log_msg('process_main: done in {} sec, #in_rec={}, #out_rec={}, saves to {}'.format(time.time() - start_time, df_main.shape, df_rstl.shape, save_fn))
    log_msg(df_rstl.describe())
    sys.stdout.flush()
    log_msg_selected_pos()

def print_text():
    now = time.strftime('%y/%m/%D-%H:%M:%S', time.localtime())
    log_msg("at %s: txt=%s" % (now, my_text))
    time.sleep(2)
    now = time.strftime('%y/%m/%D-%H:%M:%S', time.localtime())
    log_msg("at %s: txt=%s" % (now, my_text))

log_to = {}
def log_msg(msgstr, bprint=False):
    tid = threading.get_ident()
    try:
        if bprint:
            print(msgstr)
        with open(log_to[tid], 'a') as lfd:
            lfd.write('{}@{}: {} \n'.format(tid, time.strftime('%y/%m/%D-%H:%M:%S', time.localtime()), msgstr))
    except Exception as ex:
        if tid not in log_to or not os.path.exists(log_to[tid]):
            init_logger()
            with open(log_to[tid], 'a') as lfd:
                lfd.write('{}@{}: {} \n'.format(tid, time.strftime('%y/%m/%D-%H:%M:%S', time.localtime()), msgstr))
        else:
            print("log_msg: ignored exception {}".format(ex))             

def init_logger(isMain=False):
    tid = threading.get_ident()
    if tid in log_to and os.path.exists(log_to[tid]):
        print("Info: for tid {}, log file {} already created".format(tid, log_to[tid]))
    else:
        log_fn = os.path.join(LOGGER_ROOT, "rnv_{}-{}-{}.log".format('main' if isMain else 'other', proj_name, platform.node().split('.')[0][-2:])) 
        log_to[tid] = log_fn
        with open(log_to[tid], 'w') as lfd:
            lfd.write('{}@{}: Starts.. \n'.format(tid, time.strftime('%y/%m/%D-%H:%M:%S', time.localtime())))
        print("Created log file {} for tid {}".format(log_to[tid], tid))

if __name__ == '__main__':
    my_tid = threading.get_ident()
    init_logger(True)
    sys.stdout.flush()
    process_main(intpos_fn, NUM_LENGTH_LIST)
    print('Done: {}'.format(';'.join(['{} to fn {}'.format("main" if k == my_tid else "other", os.path.basename(v)) for k, v in log_to.items()])))
