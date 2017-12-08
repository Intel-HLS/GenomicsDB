import sys, os, os.path

class RepairDB(object):

    def __addLCName2run_log(self, dbpath):
        __data_handler = RunVCFData(dbpath)
        __data_handler.updateRunLogLCName()
        __data_handler.close()

    def addLCName2run_log(self, root_path):
        if os.path.isfile(root_path):
            self.__addLCName2run_log(root_path)
        else:
            dbfile = os.path.join(root_path, RunVCFData.DefaultDBName)
            if os.path.isfile(dbfile):
                self.__addLCName2run_log(dbfile)
            else:
                dirs = [ os.path.join(root_path, x) for x in os.listdir(root_path) ]
                for dn in dirs:
                    dbfile = os.path.join(dn, RunVCFData.DefaultDBName)
                    if os.path.isfile(dbfile):
                        self.__addLCName2run_log(dbfile)

    
if __name__ == '__main__':
    ws = os.path.dirname(sys.argv[0])
    sys.path.append(os.path.dirname(ws))
    from core_data import RunVCFData

    data_path = os.path.join(ws, 'vcf2tiledb-data')
    RepairDB().addLCName2run_log(data_path)
    print('DONE')