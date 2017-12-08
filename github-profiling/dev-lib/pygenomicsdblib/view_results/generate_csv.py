#! /usr/bin/python3
import sys, os, os.path
import result_handler as rh

if len(sys.argv) > 1 :
    mypath = os.path.dirname(sys.argv[0])
    result_folder = sys.argv[1] 

    print("mypath=%s, results in %s" % (mypath, result_folder))

    resultData = rh.TimeResultHandler(mypath)
    resultData.setResultPath(result_folder)

    csv_file = resultData.export2csv()

    resultData.close()
    print("DONE generated csv file %s" % csv_file)
else:
    print("Usage: python generate_csv.py result_folder")

 