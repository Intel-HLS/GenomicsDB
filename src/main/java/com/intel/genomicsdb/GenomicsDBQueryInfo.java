/*
* The MIT License (MIT)
* Copyright 2017 UCLA. License pursuant to original Intel MIT license.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy of
* this software and associated documentation files (the "Software"), to deal in
* the Software without restriction, including without limitation the rights to
* use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
* the Software, and to permit persons to whom the Software is furnished to do so,
* subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
* FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
* IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
* CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package com.intel.genomicsdb;

import java.io.IOException;
import java.io.Serializable;

/**
 * Maintain global information on query ranges which is used
 * to create GenomicsDBInputSplits
 */
class GenomicsDBQueryInfo implements Serializable {

  private long beginPosition;
  private long endPosition;

  public GenomicsDBQueryInfo(long start, long end) {
    beginPosition = start;
    endPosition = end;
  }

  public GenomicsDBQueryInfo(GenomicsDBQueryInfo copy) {
    this(copy.getBeginPosition(), copy.getEndPosition());
  }

  public long getBeginPosition() {
    return beginPosition;
  }

  public long getEndPosition() {
    return endPosition;
  }

  public boolean equals(Object obj) {
    if (obj == this) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (obj instanceof GenomicsDBQueryInfo) {
      GenomicsDBQueryInfo q = (GenomicsDBQueryInfo) obj;
      return q.getBeginPosition()==getBeginPosition() &&
              q.getEndPosition()==getEndPosition();
    }
    else {
      return false;
    }
  }
}
