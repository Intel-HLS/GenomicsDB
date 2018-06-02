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

#include "genomicsdb_jni_exception.h"
#include "genomicsdb_GenomicsDBFeatureReader.h"

#include <vector>
#include <string>

#define VERIFY_OR_THROW(X) if(!(X)) throw GenomicsDBJNIException(#X);

JNIEXPORT jboolean JNICALL Java_com_intel_genomicsdb_reader_GenomicsDBFeatureReaderJni_jniIsTileDBArray
  (JNIEnv* env, jclass currClass, jstring workspace, jstring arrayName)
{
  auto workspace_cstr = env->GetStringUTFChars(workspace, NULL);
  VERIFY_OR_THROW(workspace_cstr);
  auto array_name_cstr = env->GetStringUTFChars(arrayName, NULL);
  VERIFY_OR_THROW(array_name_cstr);

  auto return_val = array_exists(workspace_cstr, array_name_cstr);

  //Cleanup
  env->ReleaseStringUTFChars(arrayName, array_name_cstr);
  env->ReleaseStringUTFChars(workspace, workspace_cstr);

  return return_val;
}

JNIEXPORT jobjectArray JNICALL Java_com_intel_genomicsdb_reader_GenomicsDBFeatureReaderJni_jniListTileDBArrays
  (JNIEnv *env, jclass currClass, jstring workspace)
{
  auto workspace_cstr = env->GetStringUTFChars(workspace, NULL);
  VERIFY_OR_THROW(workspace_cstr);

  std::vector<std::string> array_names = get_array_names(workspace_cstr);

  jobjectArray obj_array = (jobjectArray)env->NewObjectArray(array_names.size(),
							     env->FindClass("java/lang/String"),
							     env->NewStringUTF(""));
  for(uint i=0; i<array_names.size(); i++) {
    env->SetObjectArrayElement(obj_array, i, env->NewStringUTF(array_names[i].c_str()));
  }

  env->ReleaseStringUTFChars(workspace, workspace_cstr);

  return obj_array;
}
