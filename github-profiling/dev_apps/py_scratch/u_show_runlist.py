import os
import os.path
import sys
from pprint import pprint
import pickle

proj_name = 'illmn'
root_dir = '/path/to/run_tests'

def show_runlist(all_runs):
    print("\nINFO: number of run_list file is ", len(all_runs))
    for rfn in all_runs:
        print("\nINFO: dump run_list file", rfn)
        with open(rfn, 'rb') as ifd:
            lists = pickle.load(ifd)
            for n, c in lists.items():
                print()
                print(n)
                pprint(c)

if __name__ == "__main__":
    assert len(sys.argv) > 1
    print('Start')
    flist = [os.path.join(root_dir, proj_name, 'run_list-%s' % fn) for fn in sys.argv[1:]]
    print("Dump project_name = %s, id = %d", (proj_name, sys.argv[1:]))
    show_runlist(flist)
    print('End')
