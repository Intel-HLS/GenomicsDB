#! /usr/bin/python3
# pylint: disable=abstract-method,unused-argument,no-self-use

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
class ProfilingEnv():
    def get_callsets_root(self, cust_name=None):
        print('specify the root directory that contains callsets file')
        raise NotImplementedError()

    def get_cfg_json_ws(self, cust_name=None):
        print("specify the root directory that contains json configuration files")

    def get_mapping_vid_root(self):
        print("specify the root directory that contains vid map file")
        raise NotImplementedError()

    def get_reference_root(self):
        print("specify the root directory that contains reference file")
        raise NotImplementedError()

    def get_histogram_root(self, cust_name=None):
        print("specify the root directory that contains histogram file")
        raise NotImplementedError()

    def get_interest_pos_root(self, cust_name=None):
        print("specify the root directory for interesting positions files")
        raise NotImplementedError()
    
    def get_density_pos_input(self, cust_name=None):
        print("specify the root directory")
        raise NotImplementedError()
    
    ##
    def get_hosts(self, pattern=None):
        print("specify the host list")
        raise NotImplementedError()
