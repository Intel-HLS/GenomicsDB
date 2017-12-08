import sys, os, os.path
ws = os.getcwd()
sys.path.append(ws)
from data.core_data import DataHandler
from utils.constant_def import DEFAULT_DB_NAME
from remote.get_exec_info import GenomicsExecInfo

db_path=os.path.join(ws, DEFAULT_DB_NAME)
def test_db_context():
    runid=3
    with DataHandler(db_path, False) as db:
        task_list = db.getMyRuns(runid)
    print(task_list)

def test_exec_version():
    GenomicsExecInfo(db_path).checkExecs()
    exec_fn='/home/mingrutar/cppProjects/GenomicsDB/bin/vcf2tiledb'
    version_str = GenomicsExecInfo(db_path).get_version_info(exec_fn)
    print("test_exec_version", version_str)
