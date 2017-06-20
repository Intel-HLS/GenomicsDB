/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

#ifndef JNI_MPI_INIT_H
#define JNI_MPI_INIT_H

#include "headers.h"
#ifndef DISABLE_MPI
#include <mpi.h>
#endif

class JNIMpiInit
{
  public:
    JNIMpiInit()
    {
      m_my_mpi_rank = 0;
      m_mpi_initialized = false;
    }
    void initialize()
    {
      if(!m_mpi_initialized)
      {
#ifndef DISABLE_MPI

        //Initialize MPI environment
        auto rc = MPI_Init(0, 0);
        if (rc != MPI_SUCCESS)
          printf("WARNING: MPI_Init() failed - cannot obtain MPI rank\n");
        else
        {
          //Get my world rank
          MPI_Comm_rank(MPI_COMM_WORLD, &m_my_mpi_rank);
          m_mpi_initialized = true;
        }
#endif
      }
    }
    ~JNIMpiInit()
    {
#ifndef DISABLE_MPI
      if(m_mpi_initialized)
        MPI_Finalize();
#endif
      m_mpi_initialized = false;
    }
    int get_mpi_rank() const { return m_my_mpi_rank; }
    int get_mpi_rank(const int supplied_rank) const
    {
      return (supplied_rank != 0 || !m_mpi_initialized) ? supplied_rank : m_my_mpi_rank;
    }
  private:
    int m_my_mpi_rank;
    bool m_mpi_initialized;
};

extern JNIMpiInit g_jni_mpi_init;

#endif
