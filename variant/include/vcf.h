/*  vcf.h -- VCF/BCF API functions.

    Copyright (C) 2012, 2013 Broad Institute.
    Copyright (C) 2012-2014 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */
/*
 * Stripped down version of vcf.h to re-use some values, expressions
 */
/*
    todo:
        - make the function names consistent
        - provide calls to abstract away structs as much as possible
 */

#ifndef HTSLIB_VCF_H
#define HTSLIB_VCF_H

#ifndef HTSDIR

#define BCF_VL_FIXED 0 // variable length
#define BCF_VL_VAR   1
#define BCF_VL_A     2
#define BCF_VL_G     3
#define BCF_VL_R     4


#define bcf_int32_missing    INT32_MIN
#define bcf_str_missing      0x07
extern uint32_t bcf_float_vector_end;
extern uint32_t bcf_float_missing;
typedef union
{
  uint32_t i;
  float f;
}fi_pair;
extern fi_pair bcf_float_missing_union;

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

#endif

#endif
