#! /usr/bin/python3
import platform
import datetime
import os 
import os.path
from pygenomicsdblib.utils.common import get_my_logger
from pygenomicsdblib.remote.get_exec_info import GenomicsExecInfo 

if __name__ == "__main__":
    logger = get_my_logger(os.path.basename(__file__)[:-3])
    if platform.system() != 'Windows':          # real run
        target_path = os.path.join("/", "home", "mingrutar", "cppProjects", "GenomicsDB", "bin", "vcf2tiledb")
    else:
        target_path = os.path.join("\\", "Users", "mrutarx", "bin", "Miniconda3", "Scripts", "conda.exe")
    tag = datetime.date.today().isoformat()
    exec_info_handler = GenomicsExecInfo()
    logger.info(' ==== updateBuild =====')
    exec_info_handler.updateBuild(target_path, tag)
    exec_info_handler.close()

    logger.info(' ==== checkExecs =====')
    exec_info_handler2 = GenomicsExecInfo()
    exec_info_handler2.checkExecs()

    logger.info(' ==== get_version_info =====')
    exec_info_handler3 = GenomicsExecInfo()
    version = exec_info_handler3.get_version_info(target_path)
    logger.debug(" version=%s" % version)
    exec_info_handler.close()
