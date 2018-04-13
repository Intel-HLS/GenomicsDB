#! /bin/python3

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
import re

def parse_inputs(theparser, arg_checher, required_args, optional_args, h_passthrough=None):
    theparser.add_argument('--dry-run', type=bool, nargs='?', const=True, default=False, help="dry run, generate config file only")
    def add2parser(arg, registered, isrequired):
        arg2 = re.sub(r'\W+', '', arg)
        if arg2 in registered:
            # print(arg, arg2, type(registered[arg2].the_val), registered[arg2].default_val, isrequired)
            theparser.add_argument(arg, type=type(registered[arg2].the_val), default=registered[arg2].default_val, help=registered[arg2].help_msg)
            # theparser.add_argument(arg, type=type(registered[arg2].the_val), default=registered[arg2].default_val, required=isrequired, help=registered[arg2].help_msg)
        elif isrequired:
            raise Exception('not implemented argument %s' % arg)
    if required_args:
        _ = [add2parser(arg, arg_checher.params, True) for arg in required_args]
    if optional_args:
        _ = [add2parser(arg, arg_checher.params, False) for arg in optional_args]
    if h_passthrough:
        h_passthrough(theparser)
    # my_parser.add_argument('--range', type=str, nargs='?', const='::', default='::', help="chromoseme")
    return theparser.parse_args()
