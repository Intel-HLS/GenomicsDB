#! /usr/bin/python3
''' this file used in generating callsets for illmn generator
  use 'at' to schedule next run if needed
  $ cat joba.txt
        # this is a shell file
        #clean log by  "> /var/spool/mail/mingrutar"

        hostname
        #run_cmd="source ja.bash"
        run_cmd="python3 t_check_vcfs.py"
        $run_cmd $(ls -l /data/scratch/kdatta1/illmn/*.gz | wc -l)
        ttt=$?
        echo "$run_cmd return code=$ttt"
        if [ $ttt -gt 0 ]; then
        echo "schedule next run"
        at -f joba.txt now + $ttt min             # schdule next 'at'
        at -l
        fi

$ at -f joba.txt now
'''

from subprocess import Popen

target = "TARGET_NUM"
MAX_FILE_COUNT = 3
NUM_DELTA = 5000
INITIAL_SAMPLE = 5000
MIN_INTERVAL = 1    #minimun 1 hour wait. assume 400 sample per hour

LOG_FN = "gen_callsets_copy.log"
dst_root = "/path/to/dev_data/root_dir"

def generate_hist(my_hist_fn, target_num, num_vcfs):
    print("%s: generate_hist fn=%s, #sample=%d, #available samples=%d" % \
        (datetime.datetime.now(), my_hist_fn, target_num, num_vcfs))
    proc = Popen(["python3", "/path/to/t_gen_callsets_copy.py", "-n", \
       str(target_num), "--histogram", my_hist_fn], shell=False)
    print("%s: t_gen_callsets_copy.py for %d samples, pid=%s" % (datetime.datetime.now(), target_num, proc.id))

my_num_sample = int(sys.argv[1])
my_target_num = INITIAL_SAMPLE
fn = None
next_run_int = 0
print("Start %s, init my_num_sample=%d, my_target_num=%d" % (__file__, my_num_sample, my_target_num))
hist_file_count = 0
while my_num_sample >= my_target_num and hist_file_count <= MAX_FILE_COUNT:
    hist_fn = "histogram.%s" % my_target_num
    hist_fpfn = os.path.join(dst_root, hist_fn)
    if os.path.isfile(hist_fpfn):
        print(" file %s exists, my_num_sample=%s" % (hist_fn, my_num_sample))
        my_target_num += NUM_DELTA
        hist_file_count += 1
        hist_fn = None
    else:
        print(" no file %s exists, my_num_sample=%s, call generate_hist" % (hist_fn, my_num_sample))
        generate_hist(hist_fn, my_target_num, my_num_sample)
        hist_file_count += 1
        if  hist_file_count <= MAX_FILE_COUNT:     # schedh for next
            delta_sample = my_target_num + NUM_DELTA - my_num_sample
            delta_hour = int(delta_sample/400)
            print("For next run short %d samples, check back after %d hours" % (delta_sample, delta_hour))
            if delta_hour > 0:
              next_run_int = max(MIN_INTERVAL, delta_hour)
        break
print("*** check_vcfs done next run (rc) is %d" % next_run_int)
sys.exit(next_run_int)
