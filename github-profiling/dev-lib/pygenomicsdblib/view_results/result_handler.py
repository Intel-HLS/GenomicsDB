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
import csv
from utils.constant_def import TIME_ROW_LABELS

class TestResultHandler(object):
    def __init__(self, wkspace=None ):
        self.__wspace = wkspace if wkspace else os.getcwd()
        self.__runid = None
    
    def setResultPath(self, result_path):
        self.__source = os.path.join(self.__wspace, result_path)
        sys.path.append(os.path.dirname(self.__wspace))
        import core_data
        dbpath = os.path.join(self.__source, core_data.RunVCFData.DefaultDBName)
        if os.path.isfile(dbpath):
            if hasattr(self, 'data_handler') and self.data_handler:
                self.data_handler.close() 
            self.data_handler = core_data.RunVCFData(dbpath) 
            print("INFO found genomicd loader db at %s" % dbpath )
        else:
            raise ValueError('No db @ %s' % result_path)

    def __transform4Pandas(self, inputRows) :
        ''' inputRows = [ ]
         [{k:v, k2:v2}] => [k, k2,..], [ [v, ...], [v2,...],  ]  arrya of (label, num_row)'''
        row_labels = [ x for x in inputRows[0].keys() ]
        data = []
        for i in range( len(row_labels)):
            data.append(['-'] * len(inputRows))
        for ri, row in enumerate(inputRows):
            for i, tag in enumerate(row_labels):
                data[i][ri] = row[tag]
        return row_labels, data
    
    def getRunSetting(self, runid=None):
        return self.data_handler.getRunConfigs(runid)
        
    def getRunSetting4Pandas(self, runid=None):
        ''' runid = -1 means the last run, output LOAD_x as column'''
        # confgiList [{,}] config_tag with user value or -
        my_runid, confgiList = self.data_handler.getRunConfigs(runid) 
        if confgiList:
            row_labels, data = self.__transform4Pandas(confgiList)
            col_header = [ "LOAD_%d" % (i+1) for i in range(len(confgiList))]
            pddata = [ (row_labels[i], dr) for i, dr in enumerate( data) ]
            return my_runid, pddata, col_header
        else:
            print("WARN cannot find run_id %d " % runid)

    def __get_all_result(self, runid):
        assert runid
        if runid != self.__runid or not self.__runid:
            results, cmd = self.data_handler.getAllResult(runid)
            if results:
                self.__all_results = results
                self.is_query = isinstance(results[0]['extra_info'], dict)
                self.__runid = runid
                self.genome_data = self.gexec_loader if 'vcf2tiledb' in cmd else self.gexec_query
    
    def shorten_command(self, line):
        llst = [  os.path.basename(cl) for cl in line.split() ]
        return " ".join(llst)

    def get_time_result(self, runid):
        self.__get_all_result(runid)
        col_header = []
        rtimelist = []
        for row in self.__all_results:
            row['rtime']['0'] = self.shorten_command(row['rtime']['0'])
            rtimelist.append(row['rtime'])
        row_idx, data = self.__transform4Pandas(rtimelist)
        row_labels = [ TIME_ROW_LABELS[int(i)] for  i in row_idx ]
        col_header = [ "LOAD_%d" % (i+1) for i in range(len(rtimelist))] 
        pddata = [ (row_labels[i], pd ) for i, pd in enumerate(data)]
        return pddata, col_header

    gexec_loader = {'tags' : {'fv' : 'Fetch from VCF',
                            'cc' : 'Combining Cells', 'fo' : 'Flush Output',
                            'st' : 'sections time', 'ts' : 'time in single thread phase()',
                            'tr' : 'time in read_all()'},
            'header' : ['Wall-clock time(s)', 'Cpu time(s)', 'Critical path wall-clock time(s)', 
                    'Critical path Cpu time(s)', '#critical path'] }
    gexec_query = {'tags' : { 'cf': 'GenomicsDB cell fill timer',
                                'bs' : 'bcf_t serialization',
                                'ot' : 'Operator time', 
                                'sq' : 'Sweep at query begin position', 
                                'ti' : 'TileDB iterator', 
                                'tt' : 'Total scan_and_produce_Broad_GVCF time for rank 0', 
                                'tb' : 'TileDB to buffer cell', 
                                'bc' : 'bcf_t creation time'  },
            'header' : ['Cpu time(s)', 'Wall-clock time(s)'] }

    def __get_genome_result4run(self, gendata4run):
        ll = len(gendata4run)
        row_list = [''] * ll
        for key, val in gendata4run.items():
            if key != 'op':
                if '.' in key:
                    idx = int(float(key))
                    if idx < ll:
                        row_list[idx] = val          # patch a mistake @ generator
                else:
                    row_list[int(key)] = val
        if gendata4run['op'] in self.genome_data['tags']:
            return self.genome_data['tags'][gendata4run['op']], row_list
        else:
            return gendata4run['op'], row_list

    def get_genome_results(self, runid, subidStr):
        ''' subidStr = LOAD_1,  LOAD_n '''
        self.__get_all_result(runid)
        subid = int(subidStr.split("_")[1]) - 1
        assert(subid < len(self.__all_results))
        rows = []
        for gtimes in self.__all_results[subid]['gtime'] :      # list
            rows.append(self.__get_genome_result4run(gtimes))
        return rows, self.genome_data['header']

    def get_pidstats(self, runid):
        self.__get_all_result(runid)
        pidstas = []
        for pidrow in self.__all_results:
            pidstas.append([ os.path.join(self.__source, 'stats', os.path.basename(fp)) for fp in pidrow['pidstat']])
        # pidstas = [ map(lambda fp: os.path.join(self.__source, 'stats', os.path.basename(fp)), pidrow['pidstat']) for pidrow in self.__all_results ]
        return pidstas
    
    def close(self):
        self.data_handler.close()
    
    def __csv_labels(self, confgiList):
        ldname =  next (iter (confgiList.values()))
        lc_labels = [ k for k in ldname.keys() ]
        rtime = self.__all_results[0]['rtime']      # name:val
        rtime_labels = [''] * len(rtime)
        for k,v in rtime.items():
            idx = int(k)
            rtime_labels[idx] = self.time_row_labels[idx]
        ex_labels = ["num_parallel", "qry_seg_size"] if self.is_query else ["num_parallel"]
        return ex_labels + lc_labels + rtime_labels 

    STATS_HEADS = ["%MEM", "kB_rd/s", "kB_wr/s"]
    STATS_ATTR = ["mean", "cv"]
    def write_csv_stat_labels(self, stat_writer, confgiList):
        labels = self.__csv_labels(confgiList)
        stat_labels = []
        [ stat_labels.extend([ "%s_%s" % (st, att) for att in self.STATS_ATTR ] ) for st in self.STATS_HEADS ]
        all_labels = labels + stat_labels
        stat_writer.writerow(all_labels)

    def write_csv_labels(self, csv_fd, confgiList):
        labels = self.__csv_labels(confgiList)
        gtime_labels = []
        gtimes = self.__all_results[0]['gtime']      #  30 labels num_ops x num(gtime_col_header)
        num_op = len(self.genome_data['tags'])
        perproc_count = 0
        for gt_item in gtimes:
            opname, gdata = self.__get_genome_result4run(gt_item)
            gt_op_labels = [ "%s_%s" % (opname, hname) for hname in self.genome_data['header'] ]
            gtime_labels.extend(gt_op_labels)
            perproc_count += 1
            if perproc_count % num_op == 0:
                all_labels = labels + gtime_labels
                csv_fd.write("%s\n" % ",".join(all_labels))
                csv_fd.flush()
                return

    def export2csv(self, run_dir, runid=None):
        if not runid:
            runid = self.data_handler.getLastRun()
        cmd = self.data_handler.getRunCommand(runid)
        if 'vcf2tiledb' in cmd:
            loader_run_id = runid
        else:
            loader_run_id = self.data_handler.getLoadRunId(runid)[0]
        configDict = {}
        for lcname, cfg in self.data_handler.getRunConfigsDict(loader_run_id).items():
            sorted_cfg = collections.OrderedDict( sorted(cfg.items(), key = lambda k : k[0]))   # sorted by key
            configDict[lcname] = sorted_cfg

        my_runid = runid
        assert(my_runid != None and configDict)

        print("INFO export test run %d to csv..." % my_runid)

        self.__get_all_result(my_runid)
        filename = os.path.join(self.__wspace, "csvfiles", "%s_%s-%s.csv" % (os.path.basename(run_dir), my_runid, os.path.basename(cmd)))
        csv_fd = open(filename, 'w')
        self.write_csv_labels(csv_fd, configDict)

        filename_stat = os.path.join(self.__wspace, "csvfiles", "%s_%s-%s-stat.csv" % (os.path.basename(run_dir), my_runid, os.path.basename(cmd)))
        csv_stat_fd = open(filename_stat, 'w')
        stat_writer = csv.writer(csv_stat_fd)
        self.write_csv_stat_labels(stat_writer, configDict)

        for row in self.__all_results:
            ex_data = [str(row["n_parallel"]) ]
            if self.is_query:
                ex_data.append(str(row['extra_info']['segment_size']))

            lc_data = [ v for v in configDict[row['lcname']].values() ]
            rtime = row['rtime']      # name:val
            rtime_data = [''] * len(rtime)
            for k,v in rtime.items():
                rtime_data[int(k)] = v

            for fn in row['pidstat']:
                stats_data = []
                full_fn = os.path.join(self.__source, 'stats', os.path.basename(fn))
                stat_df = pd.DataFrame.from_csv(full_fn)
                for col in self.STATS_HEADS:
                    stats_data.append(stat_df[col].mean())
                    stats_data.append(stat_df[col].std() / stats_data[-1])
                stat_row = ex_data + lc_data + rtime_data + stats_data
                stat_writer.writerow(stat_row)

            perproc_count = 0
            num_op = len(self.genome_data['tags'])
            gtime_data = []
            for gtime in row['gtime']:      # list
                opname, gdata = self.__get_genome_result4run(gtime)
                gtime_data.extend(gdata)
                perproc_count += 1
                if perproc_count % num_op == 0:
                    aRow = ex_data + lc_data + rtime_data + gtime_data
                    csv_fd.write("%s\n" % ",".join(aRow) )
                    csv_fd.flush()
                    del gtime_data[:]
        csv_fd.close()
        csv_stat_fd.close()
        return filename

def testTimeResultHandler(resultData):
    runid, confgiList = resultData.getRunSetting()

    runid, data, col = resultData.getRunSetting4Pandas()
    # test time
    time_data, col_header = resultData.get_time_result(runid)
    #test genomics data
    gen_data, col_header = resultData.get_genome_results(runid, 'LOAD_1')

    pidfiles = resultData.get_pidstats(runid)
    for pids in pidfiles:                       # each machine
        print("pifiles=%s" % ",".join(pids))

# if __name__ == '__main__':
#     mypath = os.path.dirname(sys.argv[0])
#     print("mypath=%s" % mypath)
#     resultData = TestResultHandler(mypath)

#     subdir="VDA432-compression"
#     runid = 1

#     run_dir = os.path.join("vcf2tiledb-data", subdir)
#     resultData.setResultPath(run_dir)
#     csv_file = resultData.export2csv(run_dir, runid)
#     print("csv file @ %s" % csv_file)

#     resultData.close()
#     print("DONE")
    