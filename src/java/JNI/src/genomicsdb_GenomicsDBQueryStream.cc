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

#include "genomicsdb_GenomicsDBQueryStream.h"
#include "jni_bcf_reader.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw GenomicsDBJNIException(#X);
#define GET_BCF_READER_FROM_HANDLE(X) (reinterpret_cast<JNIBCFReader*>(static_cast<std::uintptr_t>(X)))

JNIEXPORT jlong JNICALL Java_genomicsdb_GenomicsDBQueryStream_jniGenomicsDBInit
  (JNIEnv* env, jobject curr_obj, jstring loader_configuration_file, jstring query_configuration_file,
   jstring chr, jint start, jint end,
   jint rank, jlong buffer_capacity, jlong segment_size)
{
  //Java string to char*
  auto loader_configuration_file_cstr = env->GetStringUTFChars(loader_configuration_file, NULL);
  VERIFY_OR_THROW(loader_configuration_file_cstr);
  auto query_configuration_file_cstr = env->GetStringUTFChars(query_configuration_file, NULL);
  VERIFY_OR_THROW(query_configuration_file_cstr);
  auto chr_cstr = env->GetStringUTFChars(chr, NULL);
  VERIFY_OR_THROW(chr_cstr);
  //Create object
  auto output_format = "bu";
  auto bcf_reader_obj = new JNIBCFReader(loader_configuration_file_cstr, query_configuration_file_cstr,
      chr_cstr, start, end,
      rank, buffer_capacity, segment_size, output_format,
      (strcmp(output_format, "bu") == 0), false);
  //Cleanup
  env->ReleaseStringUTFChars(loader_configuration_file, loader_configuration_file_cstr);
  env->ReleaseStringUTFChars(query_configuration_file, query_configuration_file_cstr);
  env->ReleaseStringUTFChars(chr, chr_cstr);
  //Cast pointer to 64-bit int and return to Java
  return static_cast<jlong>(reinterpret_cast<std::uintptr_t>(bcf_reader_obj));
}

JNIEXPORT jlong JNICALL Java_genomicsdb_GenomicsDBQueryStream_jniGenomicsDBClose
  (JNIEnv* env, jobject curr_obj, jlong handle)
{
  auto bcf_reader_obj = GET_BCF_READER_FROM_HANDLE(handle);
  if(bcf_reader_obj) //not NULL
    delete bcf_reader_obj;
  return 0;
}

JNIEXPORT jlong JNICALL Java_genomicsdb_GenomicsDBQueryStream_jniGenomicsDBGetNumBytesAvailable
  (JNIEnv* env, jobject curr_obj, jlong handle)
{
  auto bcf_reader_obj = GET_BCF_READER_FROM_HANDLE(handle);
  return (bcf_reader_obj) ? bcf_reader_obj->get_buffer_capacity() : 0;
}

JNIEXPORT jbyte JNICALL Java_genomicsdb_GenomicsDBQueryStream_jniGenomicsDBReadNextByte
  (JNIEnv* env, jobject curr_obj, jlong handle)
{
  auto bcf_reader_obj = GET_BCF_READER_FROM_HANDLE(handle);
  return (bcf_reader_obj) ? bcf_reader_obj->read_next_byte() : -1;
}

JNIEXPORT jint JNICALL Java_genomicsdb_GenomicsDBQueryStream_jniGenomicsDBRead
  (JNIEnv* env, jobject curr_obj, jlong handle, jbyteArray java_byte_array, jint offset, jint n)
{
  auto bcf_reader_obj = GET_BCF_READER_FROM_HANDLE(handle);
  if(bcf_reader_obj == 0)
    return 0;
  auto total_num_bytes_read = 0ull;
  while(total_num_bytes_read < static_cast<uint64_t>(n) && !(bcf_reader_obj->end()))
  {
    auto& buffer_obj = bcf_reader_obj->get_read_batch();
    auto num_bytes_to_copy = std::min<size_t>(buffer_obj.get_num_remaining_bytes(), static_cast<size_t>(n)-total_num_bytes_read);
    //Handle this as a special case as read_and_advance will not advance anything if num_bytes_to_copy == 0u
    if(num_bytes_to_copy == 0u)
      num_bytes_to_copy = SIZE_MAX;     //forces jni_bcf_reader to produce the next batch of records
    else
    {
      env->SetByteArrayRegion(java_byte_array, offset+total_num_bytes_read, num_bytes_to_copy,
          reinterpret_cast<const jbyte*>(buffer_obj.get_pointer_at_read_position()));
      total_num_bytes_read += num_bytes_to_copy;
    }
    bcf_reader_obj->read_and_advance(0, 0u, num_bytes_to_copy);
  }
  return total_num_bytes_read;
}

JNIEXPORT jlong JNICALL Java_genomicsdb_GenomicsDBQueryStream_jniGenomicsDBSkip
  (JNIEnv* env, jobject curr_obj, jlong handle, jlong n)
{
  auto bcf_reader_obj = GET_BCF_READER_FROM_HANDLE(handle);
  return (bcf_reader_obj) ? bcf_reader_obj->read_and_advance(0, 0, n) : 0;
}
