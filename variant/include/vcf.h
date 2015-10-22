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
#define bcf_int32_vector_end (INT32_MIN+1)
#define bcf_str_missing      0x07
#define bcf_str_vector_end   0
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

//For bitwise comparison of float values
typedef union
{
  unsigned i;
  float f;
}fi_union;

typedef union
{
  uint64_t i;
  double d;
}di_union;

extern fi_union bcf_float_missing_union;
extern fi_union bcf_float_vector_end_union;

//Template function that return bcf missing value
template<class T>
inline T get_bcf_missing_value();

template<>
inline int get_bcf_missing_value<int>() { return bcf_int32_missing; }

template<>
inline unsigned get_bcf_missing_value<unsigned>() { return bcf_int32_missing; }

template<>
inline int64_t get_bcf_missing_value<int64_t>() { return bcf_int32_missing; }

template<>
inline uint64_t get_bcf_missing_value<uint64_t>() { return bcf_int32_missing; }

template<>
inline float get_bcf_missing_value<float>() { return bcf_float_missing_union.f; }

template<>
inline double get_bcf_missing_value<double>() { return bcf_float_missing_union.f; }

template<>
inline std::string get_bcf_missing_value<std::string>() { return ""; }

template<>
inline char get_bcf_missing_value<char>() { return bcf_str_missing; }

//Template function that return bcf vector_end value
template<class T>
inline T get_bcf_vector_end_value();

template<>
inline int get_bcf_vector_end_value<int>() { return bcf_int32_vector_end; }

template<>
inline unsigned get_bcf_vector_end_value<unsigned>() { return bcf_int32_vector_end; }

template<>
inline int64_t get_bcf_vector_end_value<int64_t>() { return bcf_int32_vector_end; }

template<>
inline uint64_t get_bcf_vector_end_value<uint64_t>() { return bcf_int32_vector_end; }

template<>
inline float get_bcf_vector_end_value<float>() { return bcf_float_vector_end_union.f; }

template<>
inline double get_bcf_vector_end_value<double>() { return bcf_float_vector_end_union.f; }

template<>
inline std::string get_bcf_vector_end_value<std::string>() { return ""; }

template<>
inline char get_bcf_vector_end_value<char>() { return bcf_str_vector_end; }

template<class T>
inline bool is_bcf_valid_value(T val) { return (val != get_bcf_missing_value<T>() && val != get_bcf_vector_end_value<T>()); }

#endif
