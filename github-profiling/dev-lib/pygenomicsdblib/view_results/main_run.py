#! /usr/bin/python3
import os.path
from result_handler import TestResultHandler

subdir="VDA432-compression"
runid = 1

if __name__ == '__main__':
    # mypath = os.path.dirname(sys.argv[0])
    # print("mypath=%s" % mypath)
    resultData = TestResultHandler()

    run_dir = os.path.join("vcf2tiledb-data", subdir)
    resultData.setResultPath(run_dir)
    csv_file = resultData.export2csv(run_dir, runid)
    print("csv file @ %s" % csv_file)

    resultData.close()
    print("DONE")
    