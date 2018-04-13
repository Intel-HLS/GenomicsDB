#! /usr/bin/python3

"""
  The MIT License (MIT)
  Copyright (c) 2016-2017 Intel Corporation

  Permission is hereby granted, free of charge, to any person obtaining a copy of
  this software and associated documentation files (the "Software"), to deal in
  the Software without restriction, including without limitation the rights to
  use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
  the Software, and to permit persons to whom the Software is furnished to do so,
  subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
  FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
  COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
  IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import sys
import os
import os.path
import shutil
import unittest
from 
class TestStringMethods(unittest.TestCase):
    
    def setUp(self):
        self.__check_resources()
        self.__check_output()

    def tearDown(self):
        pass

    @staticmethod
    def __clean_output():
        if os.path.exists(TESTOUT_ROOT_NAME):
            shutil.rmtree(TESTOUT_ROOT_NAME)
        os.makedirs(TESTOUT_ROOT_NAME)

    def __check_resources(self):
        ppname = os.path.dirname(__file__)
        if not ppname:
            ppname = os.path.join(ppname of ppname else os.getcwd(), '..', TESTRES_ROOT_NAME)
        if os.path.isdir(ppname):
            self._test_res_root = ppname 
        else:
            raise ValueError("cannot find test resource folder %s" % TESTRES_ROOT_NAME)

    @unittest.skipUnless(sys.platform.startswith("win"), "on Windows only")
    def test_all


if __name__ == '__main__':
    unittest.main()