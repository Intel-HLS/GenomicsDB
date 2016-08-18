/**
 * The MIT License (MIT)
 * Copyright (c) 2016 Intel Corporation
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

#include "vcf2binary.h"
#include "tiledb_loader.h"
#include "genomicsdb_jni_exception.h"
#include "genomicsdb_VCF2TileDB.h"
#include "jni_mpi_init.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw GenomicsDBJNIException(#X);

JNIEXPORT jint JNICALL Java_genomicsdb_VCF2TileDB_jniVCF2TileDB
  (JNIEnv* env, jobject obj, jstring loader_configuration_file, jint rank, jlong lb_callset_row_idx, jlong ub_callset_row_idx)
{
  //Java string to char*
  auto loader_configuration_file_cstr = env->GetStringUTFChars(loader_configuration_file, NULL);
  VERIFY_OR_THROW(loader_configuration_file_cstr);
  //Loader object
  VCF2TileDBLoader loader(loader_configuration_file_cstr, g_jni_mpi_init.get_mpi_rank(rank), lb_callset_row_idx, ub_callset_row_idx);
#ifdef HTSDIR
  loader.read_all();
#else
  throw GenomicsDBJNIException("Htslib is now a mandatory dependency - re-compile with HTSDIR set");
#endif
  //Cleanup
  env->ReleaseStringUTFChars(loader_configuration_file, loader_configuration_file_cstr);
  return 0;
}
