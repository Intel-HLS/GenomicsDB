#! /usr/bin/python3
import sqlite3
import os
import os.path
import datetime
import hashlib
from subprocess import check_output, CalledProcessError 
from time import sleep
from data.core_data import DataHandler
from utils.common import get_my_logger
from utils.constant_def import CHECK_VERSION_LIST, DEFAULT_DB_NAME, NO_VERSION

class GenomicsExecInfo(object):
    def __init__(self, db_path):
        self.logger = get_my_logger(self.__class__.__name__)
        self.db_path = db_path
        self.tag = None
        
    @staticmethod
    def __check_exec_hash(full_path, hashcode):
        if os.path.exists(full_path):
            checkhash = hashlib.sha256(open(full_path, 'rb').read()).hexdigest()
            return hashcode == checkhash
        else:
            raise RuntimeError("cannot find file %s" % full_path)
        
    def __insert_info(self, full_fn, need_check_version):
        hashcode = hashlib.sha256(open(full_fn, 'rb').read()).hexdigest()
        if not self.data_handler.hasVersionInfo(full_fn, hashcode):
            # we do not have it, insert it
            try:
                version = check_output([full_fn, "--version"]).decode('utf-8').strip() if need_check_version else NO_VERSION
                self.data_handler.addVersionInfo(self.tag, os.path.basename(full_fn), version, full_fn, \
                    os.path.getctime(full_fn), hashcode)
                self.logger.debug('Inserted version info of %s: %s' % (full_fn, version))
            except CalledProcessError as cpe:
                print("CalledProcessError@%s:__insert_info exception: %s" % cpe)
        else:
            self.logger.debug('Version info of %s already inserted' % full_fn)

    def __updateBuild(self, exec_path):
        if os.path.isdir(exec_path):
            for fn in next(os.walk(exec_path))[2]: 
                self.__insert_info(os.path.join(exec_path, fn), fn in CHECK_VERSION_LIST)
        else:
            self.__insert_info(exec_path, True)

    def __get_version_info(self, full_path):
        vinfo = self.data_handler.get_last_version_info(full_path)
        if vinfo and self.__check_exec_hash(full_path, vinfo[1]):
            return vinfo[0]                # we got it
        self.logger.debug("__get_version_info: no version, updated build for: %s " % full_path)
        try:
            self.__updateBuild(os.path.dirname(full_path))
        except sqlite3.OperationalError as soe:     # db is busy?
            self.logger.debug("get_version_info: " + soe)
            sleep(1)
        finally:
            version = self.data_handler.get_last_version_info(full_path)
            return version[0] if version else None
            
    def updateBuild(self, exec_path, tag=None):
        assert os.path.exists(exec_path), 'path does not exist'
        self.tag = tag if tag else datetime.date.today().isoformat()
        with DataHandler(self.db_path) as dbh:
            self.data_handler = dbh
            __updateBuild(exec_path) 
        self.data_handler = None

    def checkExecs(self):
        try:
            with DataHandler(self.db_path) as dbh:
                check_generator = dbh.checkExecsGen()
                for version, cmd, hashcode in check_generator:
                    if self.__check_exec_hash(cmd, hashcode):
                        self.logger.info("%s is valid, version='%s'" % (cmd, version))
                    else:
                        self.logger.warn("The file %s seems changed" % cmd)
        except RuntimeError as rte:
            self.logger.critical(rte)

    def get_version_info(self, full_path):
        if os.path.exists(full_path):
            with DataHandler(self.db_path) as dbh:
                self.data_handler = dbh
                version = self.__get_version_info(full_path)
            self.data_handler = None
            return version
        else:
            raise RuntimeError("cannot find file %s" % full_path)
