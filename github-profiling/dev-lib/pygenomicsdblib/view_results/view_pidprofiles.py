#! /usr/bin/python3
# plot pidstats
# generates cvs
#TODOs:
# get time results - table row= a test, 
# get genome results - table row = a stage
# duration view: bin-view x - stages, y - duration 
import sqlite3
import numpy as np
import pandas as pd
import os, os.path
import sys
import collections

if __name__ == '__main__':
    mypath = os.path.dirname(sys.argv[0])
    print("mypath=%s" % mypath)
    resultData = TimeResultHandler(mypath)

    run_dir = os.path.join("vcf2tiledb-data", "VDA413")
    resultData.setResultPath(run_dir)
    runid = 13

    # TODO: lcname, num_paral
    pidstas_files = resultData.get_pidstats(runid)
    for rowfiles in pidstas_files:
        flist = [ os.path.basename(r) for r in rowfiles]
        for 
        print('total = %d' % len(flist) )
        print(flist)


    def get_pidstats(self, runid):
        self.__get_all_result(runid)
        pidstas = []
        for pidrow in self.__all_results:
            pidstas.append([ os.path.join(self.__source, 'stats', os.path.basename(fp)) for fp in pidrow['pidstat']])
        # pidstas = [ map(lambda fp: os.path.join(self.__source, 'stats', os.path.basename(fp)), pidrow['pidstat']) for pidrow in self.__all_results ]
        return pidstas
    
    csv_file = resultData.export2csv(run_dir, runid)
    print("csv file @ %s" % csv_file)

    resultData.close()
    print("DONE")
