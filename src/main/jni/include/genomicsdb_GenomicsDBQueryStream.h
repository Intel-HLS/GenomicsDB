/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class com_intel_genomicsdb_GenomicsDBQueryStream */

#ifndef _Included_com_intel_genomicsdb_GenomicsDBQueryStream
#define _Included_com_intel_genomicsdb_GenomicsDBQueryStream
#ifdef __cplusplus
extern "C" {
#endif
#undef com_intel_genomicsdb_GenomicsDBQueryStream_MAX_SKIP_BUFFER_SIZE
#define com_intel_genomicsdb_GenomicsDBQueryStream_MAX_SKIP_BUFFER_SIZE 2048L
/*
 * Class:     com_intel_genomicsdb_GenomicsDBQueryStream
 * Method:    jniGenomicsDBInit
 * Signature: (Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;IIIJJ)J
 */
JNIEXPORT jlong JNICALL Java_com_intel_genomicsdb_GenomicsDBQueryStream_jniGenomicsDBInit
  (JNIEnv *, jobject, jstring, jstring, jstring, jint, jint, jint, jlong, jlong, jboolean, jboolean, jboolean, jboolean);

/*
 * Class:     com_intel_genomicsdb_GenomicsDBQueryStream
 * Method:    jniGenomicsDBClose
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_com_intel_genomicsdb_GenomicsDBQueryStream_jniGenomicsDBClose
  (JNIEnv *, jobject, jlong);

/*
 * Class:     com_intel_genomicsdb_GenomicsDBQueryStream
 * Method:    jniGenomicsDBGetNumBytesAvailable
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_com_intel_genomicsdb_GenomicsDBQueryStream_jniGenomicsDBGetNumBytesAvailable
  (JNIEnv *, jobject, jlong);

/*
 * Class:     com_intel_genomicsdb_GenomicsDBQueryStream
 * Method:    jniGenomicsDBReadNextByte
 * Signature: (J)B
 */
JNIEXPORT jbyte JNICALL Java_com_intel_genomicsdb_GenomicsDBQueryStream_jniGenomicsDBReadNextByte
  (JNIEnv *, jobject, jlong);

/*
 * Class:     com_intel_genomicsdb_GenomicsDBQueryStream
 * Method:    jniGenomicsDBRead
 * Signature: (J[BII)I
 */
JNIEXPORT jint JNICALL Java_com_intel_genomicsdb_GenomicsDBQueryStream_jniGenomicsDBRead
  (JNIEnv *, jobject, jlong, jbyteArray, jint, jint);

/*
 * Class:     com_intel_genomicsdb_GenomicsDBQueryStream
 * Method:    jniGenomicsDBSkip
 * Signature: (JJ)J
 */
JNIEXPORT jlong JNICALL Java_com_intel_genomicsdb_GenomicsDBQueryStream_jniGenomicsDBSkip
  (JNIEnv *, jobject, jlong, jlong);

#ifdef __cplusplus
}
#endif
#endif
