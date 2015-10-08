/*
 * Stripped down version of vcf.h to re-use some values, expressions
 * if the full htslib library is not available
 */

#ifndef HTSLIB_VCF_LITE_H
#define HTSLIB_VCF_LITE_H

#ifdef HTSDIR   //Use htslib's headers

#include "htslib/vcf.h"

#else

#define BCF_VL_FIXED 0 // variable length
#define BCF_VL_VAR   1
#define BCF_VL_A     2
#define BCF_VL_G     3
#define BCF_VL_R     4
#define BCF_VL_P     5          //ploidy


#define bcf_int32_missing    INT32_MIN
#define bcf_str_missing      0x07
extern uint32_t bcf_float_vector_end;
extern uint32_t bcf_float_missing;

static inline void bcf_float_set(float *ptr, uint32_t value)
{
    union { uint32_t i; float f; } u;
    u.i = value;
    *ptr = u.f;
}
#define bcf_float_set_vector_end(x) bcf_float_set(&(x),bcf_float_vector_end)
#define bcf_float_set_missing(x)    bcf_float_set(&(x),bcf_float_missing)
static inline int bcf_float_is_missing(float f)
{
    union { uint32_t i; float f; } u;
    u.f = f;
    return u.i==bcf_float_missing ? 1 : 0;
}
static inline int bcf_float_is_vector_end(float f)
{
    union { uint32_t i; float f; } u;
    u.f = f;
    return u.i==bcf_float_vector_end ? 1 : 0;
}

#define bcf_alleles2gt(a,b) ((a)>(b)?((a)*((a)+1)/2+(b)):((b)*((b)+1)/2+(a)))

#endif //ifdef HTSDIR

#endif
