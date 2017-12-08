#! /usr/bin/python3
# pylint: disable=broad-except
#
#
import sys
from platform import system
from astroid import MANAGER
from . run_shared import SharedTestRun

MANAGER.astroid_cache.clear()

GENERATE_CFG_ONLY = system() == 'Windows'
USE_MINI = GENERATE_CFG_ONLY

run_target = [5000, 10000, 20000, 30000, 40000]    # 5K, 10K .. 40K, 5 tiledb_ws

print("++++ Start load %s tests ++++" % __file__[2:-3])
if __name__ == '__main__':
    runner = SharedTestRun(USE_MINI)
    config_list = runner.generate_cfg(run_target)
    if config_list:
        print("... generated %d config, first two: %s" % (len(config_list), config_list[:2]))
        if not GENERATE_CFG_ONLY:
            b_cont = len(sys.argv) > 1
            runner.run_tests(config_list, b_cont)
        else:
            print("Done GENERATE_CFG_ONLY = T")
    else:
        print("Done no config file is generated")
