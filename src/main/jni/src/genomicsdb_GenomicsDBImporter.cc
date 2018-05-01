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

#include "vcf2binary.h"
#include "tiledb_loader.h"
#include "genomicsdb_importer.h"
#include "genomicsdb_jni_exception.h"
#include "genomicsdb_GenomicsDBImporter.h"
#include "jni_mpi_init.h"
#include "json_config.h"
#include "genomicsdb_vid_mapping.pb.h"
#include "genomicsdb_callsets_mapping.pb.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw GenomicsDBJNIException(#X);
#define GET_GENOMICSDB_IMPORTER_FROM_HANDLE(X) \
    (reinterpret_cast<GenomicsDBImporter*>(static_cast<std::uintptr_t>(X)))

JNIEXPORT jint JNICALL
Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniGenomicsDBImporter(
  JNIEnv* env,
  jobject obj,
  jstring loader_configuration_file,
  jint rank,
  jlong lb_callset_row_idx,
  jlong ub_callset_row_idx) {

  //Java string to char*
  auto loader_configuration_file_cstr =
      env->GetStringUTFChars(loader_configuration_file, NULL);
  VERIFY_OR_THROW(loader_configuration_file_cstr);
  //Loader object
  VCF2TileDBLoader loader(
    loader_configuration_file_cstr,
    g_jni_mpi_init.get_mpi_rank(rank),
    lb_callset_row_idx,
    ub_callset_row_idx);
#ifdef HTSDIR
  loader.read_all();
#else
  throw GenomicsDBJNIException(
    "Htslib is now a mandatory dependency - re-compile with HTSDIR set");
#endif
  //Cleanup
  env->ReleaseStringUTFChars(
    loader_configuration_file,
    loader_configuration_file_cstr);
  return 0;
}

JNIEXPORT jlong JNICALL
Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniInitializeGenomicsDBImporterObject(
  JNIEnv* env,
  jobject obj,
  jstring loader_configuration_file,
  jint rank,
  jlong lb_callset_row_idx,
  jlong ub_callset_row_idx) {

  //Java string to char*
  auto loader_configuration_file_cstr =
      env->GetStringUTFChars(loader_configuration_file, NULL);
  VERIFY_OR_THROW(loader_configuration_file_cstr);
  //GenomicsDBImporter object
  auto importer = new GenomicsDBImporter(
                        loader_configuration_file_cstr,
                        g_jni_mpi_init.get_mpi_rank(rank),
                        lb_callset_row_idx,
                        ub_callset_row_idx);
  //Cleanup
  env->ReleaseStringUTFChars(
    loader_configuration_file,
    loader_configuration_file_cstr);

  //Cast pointer to 64-bit int and return to Java
  return static_cast<jlong>(reinterpret_cast<std::uintptr_t>(importer));
}

JNIEXPORT void JNICALL
Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniAddBufferStream(
  JNIEnv* env,
  jobject obj,
  jlong handle,
  jstring stream_name,
  jboolean is_bcf,
  jlong buffer_capacity,
  jbyteArray buffer,
  jlong num_valid_bytes_in_buffer) {

  //Java string to char*
  auto stream_name_cstr = env->GetStringUTFChars(stream_name, NULL);
  VERIFY_OR_THROW(stream_name_cstr);
  jboolean is_copy = JNI_FALSE;
  //Since this function is called once per stream, performance
  //is less of a concern
  auto native_buffer_ptr = env->GetByteArrayElements(buffer, &is_copy);
  //Call importer function
  auto importer = GET_GENOMICSDB_IMPORTER_FROM_HANDLE(handle);
  importer->add_buffer_stream(
    stream_name_cstr,
    is_bcf ?
      VidFileTypeEnum::BCF_BUFFER_STREAM_TYPE :
      VidFileTypeEnum::VCF_BUFFER_STREAM_TYPE,
    buffer_capacity,
    reinterpret_cast<uint8_t*>(native_buffer_ptr),
    num_valid_bytes_in_buffer);
  //Cleanup
  env->ReleaseStringUTFChars(stream_name, stream_name_cstr);
  env->ReleaseByteArrayElements(buffer, native_buffer_ptr, 0);
}

JNIEXPORT jlong JNICALL
Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniSetupGenomicsDBLoader(
  JNIEnv* env,
  jobject obj,
  jlong handle,
  jstring buffer_stream_callset_mapping_json_string,
  jboolean usingVidMappingProtoBuf) {

  auto importer = GET_GENOMICSDB_IMPORTER_FROM_HANDLE(handle);
  //Java string to char*
  auto buffer_stream_callset_mapping_json_string_cstr =
    env->GetStringUTFChars(buffer_stream_callset_mapping_json_string, NULL);
  VERIFY_OR_THROW(buffer_stream_callset_mapping_json_string_cstr);
  importer->setup_loader(buffer_stream_callset_mapping_json_string_cstr,
      usingVidMappingProtoBuf);
  //Cleanup
  env->ReleaseStringUTFChars(
    buffer_stream_callset_mapping_json_string,
    buffer_stream_callset_mapping_json_string_cstr);
  return importer->get_max_num_buffer_stream_identifiers();
}

JNIEXPORT void JNICALL
Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniWriteDataToBufferStream(
  JNIEnv* env,
  jobject obj,
  jlong handle,
  jint buffer_stream_idx,
  jint partition_idx,
  jbyteArray buffer,
  jlong num_valid_bytes_in_buffer) {

  auto importer = GET_GENOMICSDB_IMPORTER_FROM_HANDLE(handle);
  assert(importer);
  if(importer->is_done())
    return;
  jboolean is_copy = JNI_FALSE;
  //TODO: Possible performance issues since copy is made
  auto native_buffer_ptr =
      env->GetByteArrayElements(buffer, &is_copy);
  importer->write_data_to_buffer_stream(
    buffer_stream_idx, partition_idx,
    reinterpret_cast<uint8_t*>(native_buffer_ptr),
    num_valid_bytes_in_buffer);
  //Cleanup
  env->ReleaseByteArrayElements(buffer, native_buffer_ptr, 0);
}

JNIEXPORT jboolean JNICALL
Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniImportBatch(
  JNIEnv* env,
  jobject obj,
  jlong handle,
  jlongArray exhausted_buffer_stream_identifiers) {

  auto importer = GET_GENOMICSDB_IMPORTER_FROM_HANDLE(handle);
  assert(importer);
  if(importer->is_done())
    return true;
  importer->import_batch();
  const auto& vec = importer->get_exhausted_buffer_stream_identifiers();
  //Possible buffer copy
  auto native_ebsids_ptr =
      env->GetLongArrayElements(exhausted_buffer_stream_identifiers, 0);
  for(auto i=0ull, idx=0ull;i<vec.size();++i,idx+=2u)
  {
    native_ebsids_ptr[idx] = vec[i].first;
    native_ebsids_ptr[idx+1] = vec[i].second;
  }
  //Copy number of exhausted stream identifiers in the last element
  native_ebsids_ptr[2*(importer->get_max_num_buffer_stream_identifiers())] =
      vec.size();
  //Cleanup
  env->ReleaseLongArrayElements(
    exhausted_buffer_stream_identifiers,
    native_ebsids_ptr,
    0);
  if(importer->is_done())
  {
    importer->finish();
    delete importer;
    return true;
  }
  else
    return false;
}

JNIEXPORT jstring JNICALL
Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniGetChromosomeIntervalsForColumnPartition(
  JNIEnv* env,
  jclass currClass,
  jstring loader_configuration_file,
  jint rank) {

  //Java string to char*
  auto loader_configuration_file_cstr =
      env->GetStringUTFChars(loader_configuration_file, NULL);
  VERIFY_OR_THROW(loader_configuration_file_cstr);
  auto contig_intervals =
      VidMapper::get_contig_intervals_for_column_partition(
        loader_configuration_file_cstr,
        g_jni_mpi_init.get_mpi_rank(rank),
        false);
  //Create a JSON formatted string
  rapidjson::Document json_doc;
  json_doc.SetObject();
  rapidjson::Value json_contig_intervals(rapidjson::kArrayType);
  for(const auto& curr_interval : contig_intervals)
  {
    rapidjson::Value array(rapidjson::kArrayType);
    array.PushBack(std::get<1>(curr_interval), json_doc.GetAllocator());
    array.PushBack(std::get<2>(curr_interval), json_doc.GetAllocator());
    rapidjson::Value contig_dict(rapidjson::kObjectType);
    rapidjson::Value contig_name(rapidjson::kStringType);
    contig_name.SetString(
      std::get<0>(curr_interval).c_str(),
      std::get<0>(curr_interval).length(),
      json_doc.GetAllocator());
    contig_dict.AddMember(contig_name, array, json_doc.GetAllocator());
    json_contig_intervals.PushBack(contig_dict, json_doc.GetAllocator());
  }
  json_doc.AddMember(
    "contigs",
    json_contig_intervals,
    json_doc.GetAllocator());
  rapidjson::StringBuffer buffer;
  rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
  json_doc.Accept(writer);
  //Cleanup
  env->ReleaseStringUTFChars(
    loader_configuration_file,
    loader_configuration_file_cstr);
  return env->NewStringUTF(buffer.GetString());
}


JNIEXPORT jlong JNICALL
Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniCopyVidMap
  (JNIEnv* env,
   jobject obj,
   jlong handle,
   jbyteArray vidmap_as_bytearray) {

  // Get GenomicsDBImporter object reference
  auto importer = GET_GENOMICSDB_IMPORTER_FROM_HANDLE(handle);

  VidMappingPB vidmap;
  jbyte *vidmap_elements = env->GetByteArrayElements(vidmap_as_bytearray, 0);
  int vidmap_length = env->GetArrayLength(vidmap_as_bytearray);
  try {
    vidmap.ParseFromArray(
      reinterpret_cast<jbyte*>(vidmap_elements), vidmap_length);
  } catch (std::exception e) {
    std::cerr << "Parsing error of vid maps. " <<
        "Please check VCF merged headers\n";
    throw e;
  }

  importer->setup_vidmap(&vidmap);

  //Cleanup
  env->ReleaseByteArrayElements(
        vidmap_as_bytearray,
        vidmap_elements,
        JNI_ABORT);

  //Cast pointer to 64-bit integer and return to Java
  return static_cast<jlong>(reinterpret_cast<std::uintptr_t>(importer));
} // end of jniCopyVidMap


JNIEXPORT jlong JNICALL
Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniCopyCallsetMap
  (JNIEnv* env,
   jobject obj,
   jlong handle,
   jbyteArray callsetmap_as_bytearray) {

  // Get GenomicsDBImporter object reference
  auto importer = GET_GENOMICSDB_IMPORTER_FROM_HANDLE(handle);

  CallsetMappingPB callsetmap;
  jbyte *callsetmap_elements =
      env->GetByteArrayElements(callsetmap_as_bytearray, 0);
  int callsetmap_length = env->GetArrayLength(callsetmap_as_bytearray);
  try {
    callsetmap.ParseFromArray(
        reinterpret_cast<jbyte*>(callsetmap_elements),
        callsetmap_length);
  } catch (std::exception e) {
    std::cerr << "Parsing error of callset maps. " <<
        "Please check VCF merged headers\n";
    throw e;
  }

  importer->setup_callsetmap(&callsetmap);

  //Cleanup
  env->ReleaseByteArrayElements(
        callsetmap_as_bytearray,
        callsetmap_elements,
        JNI_ABORT);

  // Cast pointer to 64-bit integer and return to Java
  return static_cast<jlong>(reinterpret_cast<std::uintptr_t>(importer));
} // end of jniCopyCallsetMap

JNIEXPORT jint JNICALL Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniCreateTileDBWorkspace
  (JNIEnv* env, jclass currClass, jstring workspace, jboolean overwrite)
{
  auto workspace_cstr = env->GetStringUTFChars(workspace, NULL);
  VERIFY_OR_THROW(workspace_cstr);
  bool should_overwrite = (overwrite == JNI_TRUE)?true:false;
  auto return_val = VCF2TileDBLoader::create_tiledb_workspace(workspace_cstr, should_overwrite);
  //Cleanup
  env->ReleaseStringUTFChars(
    workspace,
    workspace_cstr);
  return return_val;
}

JNIEXPORT jint JNICALL Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniWriteToFile
(JNIEnv* env, jclass currClass, jstring filename, jstring contents, jlong length) {
  auto filename_cstr = env->GetStringUTFChars(filename, NULL);
  VERIFY_OR_THROW(filename_cstr);
  auto contents_cstr = env->GetStringUTFChars(contents, NULL);
  VERIFY_OR_THROW(contents_cstr);
  
  auto return_val =  VCF2TileDBLoader::write_to_file(filename_cstr, (void *)contents_cstr, (size_t)length, false);

  env->ReleaseStringUTFChars(filename, filename_cstr);
  env->ReleaseStringUTFChars(contents, contents_cstr);

  return return_val;
}

JNIEXPORT jint JNICALL Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniMoveFile
(JNIEnv* env, jclass currClass, jstring source, jstring destination) {
  auto source_cstr = env->GetStringUTFChars(source, NULL);
  VERIFY_OR_THROW(source_cstr);
  auto destination_cstr = env->GetStringUTFChars(destination, NULL);
  VERIFY_OR_THROW(destination_cstr);

  auto return_val = VCF2TileDBLoader::move_file(source_cstr, destination_cstr);

  env->ReleaseStringUTFChars(source, source_cstr);
  env->ReleaseStringUTFChars(destination, destination_cstr);

  return return_val;
}

JNIEXPORT void JNICALL Java_com_intel_genomicsdb_importer_GenomicsDBImporterJni_jniConsolidateTileDBArray
  (JNIEnv* env, jclass currClass, jstring workspace, jstring array_name)
{
  auto workspace_cstr = env->GetStringUTFChars(workspace, NULL);
  VERIFY_OR_THROW(workspace_cstr);
  auto array_name_cstr = env->GetStringUTFChars(array_name, NULL);
  VERIFY_OR_THROW(array_name_cstr);
  VCF2TileDBLoader::consolidate_tiledb_array(workspace_cstr, array_name_cstr);
  //Cleanup
  env->ReleaseStringUTFChars(
      workspace,
      workspace_cstr);
  env->ReleaseStringUTFChars(
      array_name,
      array_name_cstr);
}
