import sys
import os
import os.path
import json
import logging
#from . import common
logger = logging.getLogger()

assert len(sys.argv) >= 2
callset_fn = sys.argv[1] if os.path.isfile(sys.argv[1]) else os.path.join(os.getcwd(), sys.argv[1])
assert os.path.isfile(callset_fn)

with open(callset_fn, 'r') as ifd:
    callset = json.load(ifd)
    print("number of callsets = %d" % len(callset['callsets']))
