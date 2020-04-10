/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                           License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009-2011, Willow Garage Inc., all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

/* ////////////////////////////////////////////////////////////////////
//
//  Arithmetic and logical operations: +, -, *, /, &, |, ^, ~, abs ...
//
// */

#include "precomp.hpp"

namespace libakaze
{

#if ARITHM_USE_IPP
struct IPPArithmInitializer
{
    IPPArithmInitializer(void)
    {
        ippStaticInit();
    }
};

IPPArithmInitializer ippArithmInitializer;
#endif

struct NOP {};


template<typename T, class Op, class Op8>
void vBinOp8(const T* src1, size_t step1, const T* src2, size_t step2, T* dst, size_t step, Size sz)
{
#if CV_SSE2
    Op8 op8;
#endif
    Op op;

    for( ; sz.height--; src1 += step1/sizeof(src1[0]),
                        src2 += step2/sizeof(src2[0]),
                        dst += step/sizeof(dst[0]) )
    {
        int x = 0;

    #if CV_SSE2
        if( USE_SSE2 )
        {
            for( ; x <= sz.width - 32; x += 32 )
            {
                __m128i r0 = _mm_loadu_si128((const __m128i*)(src1 + x));
                __m128i r1 = _mm_loadu_si128((const __m128i*)(src1 + x + 16));
                r0 = op8(r0,_mm_loadu_si128((const __m128i*)(src2 + x)));
                r1 = op8(r1,_mm_loadu_si128((const __m128i*)(src2 + x + 16)));
                _mm_storeu_si128((__m128i*)(dst + x), r0);
                _mm_storeu_si128((__m128i*)(dst + x + 16), r1);
            }
            for( ; x <= sz.width - 8; x += 8 )
            {
                __m128i r0 = _mm_loadl_epi64((const __m128i*)(src1 + x));
                r0 = op8(r0,_mm_loadl_epi64((const __m128i*)(src2 + x)));
                _mm_storel_epi64((__m128i*)(dst + x), r0);
            }
        }
    #endif
#if CV_ENABLE_UNROLLED
        for( ; x <= sz.width - 4; x += 4 )
        {
            T v0 = op(src1[x], src2[x]);
            T v1 = op(src1[x+1], src2[x+1]);
            dst[x] = v0; dst[x+1] = v1;
            v0 = op(src1[x+2], src2[x+2]);
            v1 = op(src1[x+3], src2[x+3]);
            dst[x+2] = v0; dst[x+3] = v1;
        }
#endif
        for( ; x < sz.width; x++ )
            dst[x] = op(src1[x], src2[x]);
    }
}

template<typename T, class Op, class Op16>
void vBinOp16(const T* src1, size_t step1, const T* src2, size_t step2,
              T* dst, size_t step, Size sz)
{
#if CV_SSE2
    Op16 op16;
#endif
    Op op;

    for( ; sz.height--; src1 += step1/sizeof(src1[0]),
        src2 += step2/sizeof(src2[0]),
        dst += step/sizeof(dst[0]) )
    {
        int x = 0;

    #if CV_SSE2
        if( USE_SSE2 )
        {
            for( ; x <= sz.width - 16; x += 16 )
            {
                __m128i r0 = _mm_loadu_si128((const __m128i*)(src1 + x));
                __m128i r1 = _mm_loadu_si128((const __m128i*)(src1 + x + 8));
                r0 = op16(r0,_mm_loadu_si128((const __m128i*)(src2 + x)));
                r1 = op16(r1,_mm_loadu_si128((const __m128i*)(src2 + x + 8)));
                _mm_storeu_si128((__m128i*)(dst + x), r0);
                _mm_storeu_si128((__m128i*)(dst + x + 8), r1);
            }
            for( ; x <= sz.width - 4; x += 4 )
            {
                __m128i r0 = _mm_loadl_epi64((const __m128i*)(src1 + x));
                r0 = op16(r0,_mm_loadl_epi64((const __m128i*)(src2 + x)));
                _mm_storel_epi64((__m128i*)(dst + x), r0);
            }
        }
        else
    #endif

        for( ; x <= sz.width - 4; x += 4 )
        {
            T v0 = op(src1[x], src2[x]);
            T v1 = op(src1[x+1], src2[x+1]);
            dst[x] = v0; dst[x+1] = v1;
            v0 = op(src1[x+2], src2[x+2]);
            v1 = op(src1[x+3], src2[x+3]);
            dst[x+2] = v0; dst[x+3] = v1;
        }

        for( ; x < sz.width; x++ )
            dst[x] = op(src1[x], src2[x]);
    }
}


template<class Op, class Op32>
void vBinOp32s(const int* src1, size_t step1, const int* src2, size_t step2,
               int* dst, size_t step, Size sz)
{
#if CV_SSE2
    Op32 op32;
#endif
    Op op;

    for( ; sz.height--; src1 += step1/sizeof(src1[0]),
        src2 += step2/sizeof(src2[0]),
        dst += step/sizeof(dst[0]) )
    {
        int x = 0;

#if CV_SSE2
        if( USE_SSE2 )
        {
            if( (((size_t)src1|(size_t)src2|(size_t)dst)&15) == 0 )
                for( ; x <= sz.width - 8; x += 8 )
                {
                    __m128i r0 = _mm_load_si128((const __m128i*)(src1 + x));
                    __m128i r1 = _mm_load_si128((const __m128i*)(src1 + x + 4));
                    r0 = op32(r0,_mm_load_si128((const __m128i*)(src2 + x)));
                    r1 = op32(r1,_mm_load_si128((const __m128i*)(src2 + x + 4)));
                    _mm_store_si128((__m128i*)(dst + x), r0);
                    _mm_store_si128((__m128i*)(dst + x + 4), r1);
                }
            else
                for( ; x <= sz.width - 8; x += 8 )
                {
                    __m128i r0 = _mm_loadu_si128((const __m128i*)(src1 + x));
                    __m128i r1 = _mm_loadu_si128((const __m128i*)(src1 + x + 4));
                    r0 = op32(r0,_mm_loadu_si128((const __m128i*)(src2 + x)));
                    r1 = op32(r1,_mm_loadu_si128((const __m128i*)(src2 + x + 4)));
                    _mm_storeu_si128((__m128i*)(dst + x), r0);
                    _mm_storeu_si128((__m128i*)(dst + x + 4), r1);
                }
        }
#endif
#if CV_ENABLE_UNROLLED
        for( ; x <= sz.width - 4; x += 4 )
        {
            int v0 = op(src1[x], src2[x]);
            int v1 = op(src1[x+1], src2[x+1]);
            dst[x] = v0; dst[x+1] = v1;
            v0 = op(src1[x+2], src2[x+2]);
            v1 = op(src1[x+3], src2[x+3]);
            dst[x+2] = v0; dst[x+3] = v1;
        }
#endif
        for( ; x < sz.width; x++ )
            dst[x] = op(src1[x], src2[x]);
    }
}


template<class Op, class Op32>
void vBinOp32f(const float* src1, size_t step1, const float* src2, size_t step2,
               float* dst, size_t step, Size sz)
{
#if CV_SSE2
    Op32 op32;
#endif
    Op op;

    for( ; sz.height--; src1 += step1/sizeof(src1[0]),
        src2 += step2/sizeof(src2[0]),
        dst += step/sizeof(dst[0]) )
    {
        int x = 0;

    #if CV_SSE2
        if( USE_SSE2 )
        {
            if( (((size_t)src1|(size_t)src2|(size_t)dst)&15) == 0 )
                for( ; x <= sz.width - 8; x += 8 )
                {
                    __m128 r0 = _mm_load_ps(src1 + x);
                    __m128 r1 = _mm_load_ps(src1 + x + 4);
                    r0 = op32(r0,_mm_load_ps(src2 + x));
                    r1 = op32(r1,_mm_load_ps(src2 + x + 4));
                    _mm_store_ps(dst + x, r0);
                    _mm_store_ps(dst + x + 4, r1);
                }
            else
                for( ; x <= sz.width - 8; x += 8 )
                {
                    __m128 r0 = _mm_loadu_ps(src1 + x);
                    __m128 r1 = _mm_loadu_ps(src1 + x + 4);
                    r0 = op32(r0,_mm_loadu_ps(src2 + x));
                    r1 = op32(r1,_mm_loadu_ps(src2 + x + 4));
                    _mm_storeu_ps(dst + x, r0);
                    _mm_storeu_ps(dst + x + 4, r1);
                }
        }
    #endif
#if CV_ENABLE_UNROLLED
        for( ; x <= sz.width - 4; x += 4 )
        {
            float v0 = op(src1[x], src2[x]);
            float v1 = op(src1[x+1], src2[x+1]);
            dst[x] = v0; dst[x+1] = v1;
            v0 = op(src1[x+2], src2[x+2]);
            v1 = op(src1[x+3], src2[x+3]);
            dst[x+2] = v0; dst[x+3] = v1;
        }
#endif
        for( ; x < sz.width; x++ )
            dst[x] = op(src1[x], src2[x]);
    }
}

template<class Op, class Op64>
void vBinOp64f(const double* src1, size_t step1, const double* src2, size_t step2,
               double* dst, size_t step, Size sz)
{
#if CV_SSE2
    Op64 op64;
#endif
    Op op;

    for( ; sz.height--; src1 += step1/sizeof(src1[0]),
        src2 += step2/sizeof(src2[0]),
        dst += step/sizeof(dst[0]) )
    {
        int x = 0;

    #if CV_SSE2
        if( USE_SSE2 && (((size_t)src1|(size_t)src2|(size_t)dst)&15) == 0 )
            for( ; x <= sz.width - 4; x += 4 )
            {
                __m128d r0 = _mm_load_pd(src1 + x);
                __m128d r1 = _mm_load_pd(src1 + x + 2);
                r0 = op64(r0,_mm_load_pd(src2 + x));
                r1 = op64(r1,_mm_load_pd(src2 + x + 2));
                _mm_store_pd(dst + x, r0);
                _mm_store_pd(dst + x + 2, r1);
            }
        else
    #endif
        for( ; x <= sz.width - 4; x += 4 )
        {
            double v0 = op(src1[x], src2[x]);
            double v1 = op(src1[x+1], src2[x+1]);
            dst[x] = v0; dst[x+1] = v1;
            v0 = op(src1[x+2], src2[x+2]);
            v1 = op(src1[x+3], src2[x+3]);
            dst[x+2] = v0; dst[x+3] = v1;
        }

        for( ; x < sz.width; x++ )
            dst[x] = op(src1[x], src2[x]);
    }
}

#if CV_SSE2

struct _VAdd8u { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_adds_epu8(a,b); }};
struct _VSub8u { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_subs_epu8(a,b); }};
struct _VMin8u { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_min_epu8(a,b); }};
struct _VMax8u { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_max_epu8(a,b); }};
struct _VAbsDiff8u
{
    __m128i operator()(const __m128i& a, const __m128i& b) const
    { return _mm_add_epi8(_mm_subs_epu8(a,b),_mm_subs_epu8(b,a)); }
};

struct _VAdd8s { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_adds_epi8(a,b); }};
struct _VSub8s { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_subs_epi8(a,b); }};
struct _VMin8s
{
    __m128i operator()(const __m128i& a, const __m128i& b) const
    {
        __m128i m = _mm_cmpgt_epi8(a, b);
        return _mm_xor_si128(a, _mm_and_si128(_mm_xor_si128(a, b), m));
    }
};
struct _VMax8s
{
    __m128i operator()(const __m128i& a, const __m128i& b) const
    {
        __m128i m = _mm_cmpgt_epi8(b, a);
        return _mm_xor_si128(a, _mm_and_si128(_mm_xor_si128(a, b), m));
    }
};
struct _VAbsDiff8s
{
    __m128i operator()(const __m128i& a, const __m128i& b) const
    {
        __m128i d = _mm_subs_epi8(a, b);
        __m128i m = _mm_cmpgt_epi8(b, a);
        return _mm_subs_epi8(_mm_xor_si128(d, m), m);
    }
};

struct _VAdd16u { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_adds_epu16(a,b); }};
struct _VSub16u { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_subs_epu16(a,b); }};
struct _VMin16u
{
    __m128i operator()(const __m128i& a, const __m128i& b) const
    { return _mm_subs_epu16(a,_mm_subs_epu16(a,b)); }
};
struct _VMax16u
{
    __m128i operator()(const __m128i& a, const __m128i& b) const
    { return _mm_adds_epu16(_mm_subs_epu16(a,b),b); }
};
struct _VAbsDiff16u
{
    __m128i operator()(const __m128i& a, const __m128i& b) const
    { return _mm_add_epi16(_mm_subs_epu16(a,b),_mm_subs_epu16(b,a)); }
};

struct _VAdd16s { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_adds_epi16(a,b); }};
struct _VSub16s { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_subs_epi16(a,b); }};
struct _VMin16s { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_min_epi16(a,b); }};
struct _VMax16s { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_max_epi16(a,b); }};
struct _VAbsDiff16s
{
    __m128i operator()(const __m128i& a, const __m128i& b) const
    {
        __m128i M = _mm_max_epi16(a,b), m = _mm_min_epi16(a,b);
        return _mm_subs_epi16(M, m);
    }
};

struct _VAdd32s { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_add_epi32(a,b); }};
struct _VSub32s { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_sub_epi32(a,b); }};
struct _VMin32s
{
    __m128i operator()(const __m128i& a, const __m128i& b) const
    {
        __m128i m = _mm_cmpgt_epi32(a, b);
        return _mm_xor_si128(a, _mm_and_si128(_mm_xor_si128(a, b), m));
    }
};
struct _VMax32s
{
    __m128i operator()(const __m128i& a, const __m128i& b) const
    {
        __m128i m = _mm_cmpgt_epi32(b, a);
        return _mm_xor_si128(a, _mm_and_si128(_mm_xor_si128(a, b), m));
    }
};
struct _VAbsDiff32s
{
    __m128i operator()(const __m128i& a, const __m128i& b) const
    {
        __m128i d = _mm_sub_epi32(a, b);
        __m128i m = _mm_cmpgt_epi32(b, a);
        return _mm_sub_epi32(_mm_xor_si128(d, m), m);
    }
};

struct _VAdd32f { __m128 operator()(const __m128& a, const __m128& b) const { return _mm_add_ps(a,b); }};
struct _VSub32f { __m128 operator()(const __m128& a, const __m128& b) const { return _mm_sub_ps(a,b); }};
struct _VMin32f { __m128 operator()(const __m128& a, const __m128& b) const { return _mm_min_ps(a,b); }};
struct _VMax32f { __m128 operator()(const __m128& a, const __m128& b) const { return _mm_max_ps(a,b); }};
static int CV_DECL_ALIGNED(16) v32f_absmask[] = { 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff };
struct _VAbsDiff32f
{
    __m128 operator()(const __m128& a, const __m128& b) const
    {
        return _mm_and_ps(_mm_sub_ps(a,b), *(const __m128*)v32f_absmask);
    }
};

struct _VAdd64f { __m128d operator()(const __m128d& a, const __m128d& b) const { return _mm_add_pd(a,b); }};
struct _VSub64f { __m128d operator()(const __m128d& a, const __m128d& b) const { return _mm_sub_pd(a,b); }};
struct _VMin64f { __m128d operator()(const __m128d& a, const __m128d& b) const { return _mm_min_pd(a,b); }};
struct _VMax64f { __m128d operator()(const __m128d& a, const __m128d& b) const { return _mm_max_pd(a,b); }};

static unsigned int CV_DECL_ALIGNED(16) v64f_absmask[] = { 0xffffffff, 0x7fffffff, 0xffffffff, 0x7fffffff };
struct _VAbsDiff64f
{
    __m128d operator()(const __m128d& a, const __m128d& b) const
    {
        return _mm_and_pd(_mm_sub_pd(a,b), *(const __m128d*)v64f_absmask);
    }
};

struct _VAnd8u { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_and_si128(a,b); }};
struct _VOr8u  { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_or_si128(a,b); }};
struct _VXor8u { __m128i operator()(const __m128i& a, const __m128i& b) const { return _mm_xor_si128(a,b); }};
struct _VNot8u { __m128i operator()(const __m128i& a, const __m128i&) const { return _mm_xor_si128(_mm_set1_epi32(-1),a); }};

#endif

#if CV_SSE2
#define IF_SIMD(op) op
#else
#define IF_SIMD(op) NOP
#endif

template<> inline uchar OpAdd<uchar>::operator ()(uchar a, uchar b) const
{ return CV_FAST_CAST_8U(a + b); }
template<> inline uchar OpSub<uchar>::operator ()(uchar a, uchar b) const
{ return CV_FAST_CAST_8U(a - b); }

template<typename T> struct OpAbsDiff
{
    typedef T type1;
    typedef T type2;
    typedef T rtype;
    T operator()(T a, T b) const { return (T)std::abs(a - b); }
};

template<> inline short OpAbsDiff<short>::operator ()(short a, short b) const
{ return saturate_cast<short>(std::abs(a - b)); }

template<> inline schar OpAbsDiff<schar>::operator ()(schar a, schar b) const
{ return saturate_cast<schar>(std::abs(a - b)); }

template<typename T, typename WT=T> struct OpAbsDiffS
{
    typedef T type1;
    typedef WT type2;
    typedef T rtype;
    T operator()(T a, WT b) const { return saturate_cast<T>(std::abs(a - b)); }
};

template<typename T> struct OpAnd
{
    typedef T type1;
    typedef T type2;
    typedef T rtype;
    T operator()( T a, T b ) const { return a & b; }
};

template<typename T> struct OpOr
{
    typedef T type1;
    typedef T type2;
    typedef T rtype;
    T operator()( T a, T b ) const { return a | b; }
};

template<typename T> struct OpXor
{
    typedef T type1;
    typedef T type2;
    typedef T rtype;
    T operator()( T a, T b ) const { return a ^ b; }
};

template<typename T> struct OpNot
{
    typedef T type1;
    typedef T type2;
    typedef T rtype;
    T operator()( T a, T ) const { return ~a; }
};

static inline void fixSteps(Size sz, size_t elemSize, size_t& step1, size_t& step2, size_t& step)
{
    if( sz.height == 1 )
        step1 = step2 = step = sz.width*elemSize;
}

static void add8u( const uchar* src1, size_t step1,
                   const uchar* src2, size_t step2,
                   uchar* dst, size_t step, Size sz, void* )
{
    IF_IPP(fixSteps(sz, sizeof(dst[0]), step1, step2, step);
           ippiAdd_8u_C1RSfs(src1, (int)step1, src2, (int)step2, dst, (int)step, (IppiSize&)sz, 0),
           (vBinOp8<uchar, OpAdd<uchar>, IF_SIMD(_VAdd8u)>(src1, step1, src2, step2, dst, step, sz)));
}

static void add8s( const schar* src1, size_t step1,
                   const schar* src2, size_t step2,
                   schar* dst, size_t step, Size sz, void* )
{
    vBinOp8<schar, OpAdd<schar>, IF_SIMD(_VAdd8s)>(src1, step1, src2, step2, dst, step, sz);
}

static void add16u( const ushort* src1, size_t step1,
                    const ushort* src2, size_t step2,
                    ushort* dst, size_t step, Size sz, void* )
{
    IF_IPP(fixSteps(sz, sizeof(dst[0]), step1, step2, step);
           ippiAdd_16u_C1RSfs(src1, (int)step1, src2, (int)step2, dst, (int)step, (IppiSize&)sz, 0),
            (vBinOp16<ushort, OpAdd<ushort>, IF_SIMD(_VAdd16u)>(src1, step1, src2, step2, dst, step, sz)));
}

static void add16s( const short* src1, size_t step1,
                    const short* src2, size_t step2,
                    short* dst, size_t step, Size sz, void* )
{
    IF_IPP(fixSteps(sz, sizeof(dst[0]), step1, step2, step);
           ippiAdd_16s_C1RSfs(src1, (int)step1, src2, (int)step2, dst, (int)step, (IppiSize&)sz, 0),
           (vBinOp16<short, OpAdd<short>, IF_SIMD(_VAdd16s)>(src1, step1, src2, step2, dst, step, sz)));
}

static void add32s( const int* src1, size_t step1,
                    const int* src2, size_t step2,
                    int* dst, size_t step, Size sz, void* )
{
    vBinOp32s<OpAdd<int>, IF_SIMD(_VAdd32s)>(src1, step1, src2, step2, dst, step, sz);
}

static void add32f( const float* src1, size_t step1,
                    const float* src2, size_t step2,
                    float* dst, size_t step, Size sz, void* )
{
    IF_IPP(fixSteps(sz, sizeof(dst[0]), step1, step2, step);
           ippiAdd_32f_C1R(src1, (int)step1, src2, (int)step2, dst, (int)step, (IppiSize&)sz),
           (vBinOp32f<OpAdd<float>, IF_SIMD(_VAdd32f)>(src1, step1, src2, step2, dst, step, sz)));
}

static void add64f( const double* src1, size_t step1,
                    const double* src2, size_t step2,
                    double* dst, size_t step, Size sz, void* )
{
    vBinOp64f<OpAdd<double>, IF_SIMD(_VAdd64f)>(src1, step1, src2, step2, dst, step, sz);
}

static void sub8u( const uchar* src1, size_t step1,
                   const uchar* src2, size_t step2,
                   uchar* dst, size_t step, Size sz, void* )
{
    IF_IPP(fixSteps(sz, sizeof(dst[0]), step1, step2, step);
           ippiSub_8u_C1RSfs(src2, (int)step2, src1, (int)step1, dst, (int)step, (IppiSize&)sz, 0),
           (vBinOp8<uchar, OpSub<uchar>, IF_SIMD(_VSub8u)>(src1, step1, src2, step2, dst, step, sz)));
}

static void sub8s( const schar* src1, size_t step1,
                   const schar* src2, size_t step2,
                   schar* dst, size_t step, Size sz, void* )
{
    vBinOp8<schar, OpSub<schar>, IF_SIMD(_VSub8s)>(src1, step1, src2, step2, dst, step, sz);
}

static void sub16u( const ushort* src1, size_t step1,
                    const ushort* src2, size_t step2,
                    ushort* dst, size_t step, Size sz, void* )
{
    IF_IPP(fixSteps(sz, sizeof(dst[0]), step1, step2, step);
           ippiSub_16u_C1RSfs(src2, (int)step2, src1, (int)step1, dst, (int)step, (IppiSize&)sz, 0),
           (vBinOp16<ushort, OpSub<ushort>, IF_SIMD(_VSub16u)>(src1, step1, src2, step2, dst, step, sz)));
}

static void sub16s( const short* src1, size_t step1,
                    const short* src2, size_t step2,
                    short* dst, size_t step, Size sz, void* )
{
    IF_IPP(fixSteps(sz, sizeof(dst[0]), step1, step2, step);
           ippiSub_16s_C1RSfs(src2, (int)step2, src1, (int)step1, dst, (int)step, (IppiSize&)sz, 0),
           (vBinOp16<short, OpSub<short>, IF_SIMD(_VSub16s)>(src1, step1, src2, step2, dst, step, sz)));
}

static void sub32s( const int* src1, size_t step1,
                    const int* src2, size_t step2,
                    int* dst, size_t step, Size sz, void* )
{
    vBinOp32s<OpSub<int>, IF_SIMD(_VSub32s)>(src1, step1, src2, step2, dst, step, sz);
}

static void sub32f( const float* src1, size_t step1,
                   const float* src2, size_t step2,
                   float* dst, size_t step, Size sz, void* )
{
    IF_IPP(fixSteps(sz, sizeof(dst[0]), step1, step2, step);
           ippiSub_32f_C1R(src2, (int)step2, src1, (int)step1, dst, (int)step, (IppiSize&)sz),
           (vBinOp32f<OpSub<float>, IF_SIMD(_VSub32f)>(src1, step1, src2, step2, dst, step, sz)));
}

static void sub64f( const double* src1, size_t step1,
                    const double* src2, size_t step2,
                    double* dst, size_t step, Size sz, void* )
{
    vBinOp64f<OpSub<double>, IF_SIMD(_VSub64f)>(src1, step1, src2, step2, dst, step, sz);
}






static void and8u( const uchar* src1, size_t step1,
                   const uchar* src2, size_t step2,
                   uchar* dst, size_t step, Size sz, void* )
{
    IF_IPP(fixSteps(sz, sizeof(dst[0]), step1, step2, step);
           ippiAnd_8u_C1R(src1, (int)step1, src2, (int)step2, dst, (int)step, (IppiSize&)sz),
           (vBinOp8<uchar, OpAnd<uchar>, IF_SIMD(_VAnd8u)>(src1, step1, src2, step2, dst, step, sz)));
}

static void or8u( const uchar* src1, size_t step1,
                  const uchar* src2, size_t step2,
                  uchar* dst, size_t step, Size sz, void* )
{
    IF_IPP(fixSteps(sz, sizeof(dst[0]), step1, step2, step);
           ippiOr_8u_C1R(src1, (int)step1, src2, (int)step2, dst, (int)step, (IppiSize&)sz),
           (vBinOp8<uchar, OpOr<uchar>, IF_SIMD(_VOr8u)>(src1, step1, src2, step2, dst, step, sz)));
}

static void xor8u( const uchar* src1, size_t step1,
                   const uchar* src2, size_t step2,
                   uchar* dst, size_t step, Size sz, void* )
{
    IF_IPP(fixSteps(sz, sizeof(dst[0]), step1, step2, step);
           ippiXor_8u_C1R(src1, (int)step1, src2, (int)step2, dst, (int)step, (IppiSize&)sz),
           (vBinOp8<uchar, OpXor<uchar>, IF_SIMD(_VXor8u)>(src1, step1, src2, step2, dst, step, sz)));
}


/****************************************************************************************\
*                                   logical operations                                   *
\****************************************************************************************/

void convertAndUnrollScalar( const Mat& sc, int buftype, uchar* scbuf, size_t blocksize )
{
    int scn = (int)sc.total(), cn = CV_MAT_CN(buftype);
    size_t esz = CV_ELEM_SIZE(buftype);
    getConvertFunc(sc.depth(), buftype)(sc.data, 0, 0, 0, scbuf, 0, Size(std::min(cn, scn), 1), 0);
    // unroll the scalar
    if( scn < cn )
    {
        CV_Assert( scn == 1 );
        size_t esz1 = CV_ELEM_SIZE1(buftype);
        for( size_t i = esz1; i < esz; i++ )
            scbuf[i] = scbuf[i - esz1];
    }
    for( size_t i = esz; i < blocksize*esz; i++ )
        scbuf[i] = scbuf[i - esz];
}

static void binary_op(InputArray _src1, InputArray _src2, OutputArray _dst,
               InputArray _mask, const BinaryFunc* tab, bool bitwise)
{
    int kind1 = _src1.kind(), kind2 = _src2.kind();
    Mat src1 = _src1.getMat(), src2 = _src2.getMat();
    bool haveMask = !_mask.empty(), haveScalar = false;
    BinaryFunc func;
    int c;

    if( src1.dims <= 2 && src2.dims <= 2 && kind1 == kind2 &&
        src1.size() == src2.size() && src1.type() == src2.type() && !haveMask )
    {
        _dst.create(src1.size(), src1.type());
        Mat dst = _dst.getMat();
        if( bitwise )
        {
            func = *tab;
            c = (int)src1.elemSize();
        }
        else
        {
            func = tab[src1.depth()];
            c = src1.channels();
        }

        Size sz = getContinuousSize(src1, src2, dst);
        size_t len = sz.width*(size_t)c;
        if( len == (size_t)(int)len )
        {
            sz.width = (int)len;
            func(src1.data, src1.step, src2.data, src2.step, dst.data, dst.step, sz, 0);
            return;
        }
    }

    if( (kind1 == _InputArray::MATX) + (kind2 == _InputArray::MATX) == 1 ||
        src1.size != src2.size || src1.type() != src2.type() )
    {
        if( checkScalar(src1, src2.type(), kind1, kind2) )
            // src1 is a scalar; swap it with src2
            swap(src1, src2);
        else if( !checkScalar(src2, src1.type(), kind2, kind1) )
            CV_Error( CV_StsUnmatchedSizes,
                      "The operation is neither 'array op array' (where arrays have the same size and type), "
                      "nor 'array op scalar', nor 'scalar op array'" );
        haveScalar = true;
    }

    size_t esz = src1.elemSize();
    size_t blocksize0 = (BLOCK_SIZE + esz-1)/esz;
    int cn = src1.channels();
    BinaryFunc copymask = 0;
    Mat mask;
    bool reallocate = false;

    if( haveMask )
    {
        mask = _mask.getMat();
        CV_Assert( (mask.type() == CV_8UC1 || mask.type() == CV_8SC1) );
        CV_Assert( mask.size == src1.size );
        copymask = getCopyMaskFunc(esz);
        Mat tdst = _dst.getMat();
        reallocate = tdst.size != src1.size || tdst.type() != src1.type();
    }

    AutoBuffer<uchar> _buf;
    uchar *scbuf = 0, *maskbuf = 0;

    _dst.create(src1.dims, src1.size, src1.type());
    Mat dst = _dst.getMat();

    // if this is mask operation and dst has been reallocated,
    // we have to
    if( haveMask && reallocate )
        dst = Scalar::all(0);

    if( bitwise )
    {
        func = *tab;
        c = (int)esz;
    }
    else
    {
        func = tab[src1.depth()];
        c = cn;
    }

    if( !haveScalar )
    {
        const Mat* arrays[] = { &src1, &src2, &dst, &mask, 0 };
        uchar* ptrs[4];

        NAryMatIterator it(arrays, ptrs);
        size_t total = it.size, blocksize = total;

        if( blocksize*c > INT_MAX )
            blocksize = INT_MAX/c;

        if( haveMask )
        {
            blocksize = std::min(blocksize, blocksize0);
            _buf.allocate(blocksize*esz);
            maskbuf = _buf;
        }

        for( size_t i = 0; i < it.nplanes; i++, ++it )
        {
            for( size_t j = 0; j < total; j += blocksize )
            {
                int bsz = (int)MIN(total - j, blocksize);

                func( ptrs[0], 0, ptrs[1], 0, haveMask ? maskbuf : ptrs[2], 0, Size(bsz*c, 1), 0 );
                if( haveMask )
                {
                    copymask( maskbuf, 0, ptrs[3], 0, ptrs[2], 0, Size(bsz, 1), &esz );
                    ptrs[3] += bsz;
                }

                bsz *= (int)esz;
                ptrs[0] += bsz; ptrs[1] += bsz; ptrs[2] += bsz;
            }
        }
    }
    else
    {
        const Mat* arrays[] = { &src1, &dst, &mask, 0 };
        uchar* ptrs[3];

        NAryMatIterator it(arrays, ptrs);
        size_t total = it.size, blocksize = std::min(total, blocksize0);

        _buf.allocate(blocksize*(haveMask ? 2 : 1)*esz + 32);
        scbuf = _buf;
        maskbuf = alignPtr(scbuf + blocksize*esz, 16);

        convertAndUnrollScalar( src2, src1.type(), scbuf, blocksize);

        for( size_t i = 0; i < it.nplanes; i++, ++it )
        {
            for( size_t j = 0; j < total; j += blocksize )
            {
                int bsz = (int)MIN(total - j, blocksize);

                func( ptrs[0], 0, scbuf, 0, haveMask ? maskbuf : ptrs[1], 0, Size(bsz*c, 1), 0 );
                if( haveMask )
                {
                    copymask( maskbuf, 0, ptrs[2], 0, ptrs[1], 0, Size(bsz, 1), &esz );
                    ptrs[2] += bsz;
                }

                bsz *= (int)esz;
                ptrs[0] += bsz; ptrs[1] += bsz;
            }
        }
    }
}


}

void libakaze::bitwise_and(InputArray a, InputArray b, OutputArray c, InputArray mask)
{
    BinaryFunc f = (BinaryFunc)GET_OPTIMIZED(and8u);
    binary_op(a, b, c, mask, &f, true);
}

void libakaze::bitwise_or(InputArray a, InputArray b, OutputArray c, InputArray mask)
{
    BinaryFunc f = (BinaryFunc)GET_OPTIMIZED(or8u);
    binary_op(a, b, c, mask, &f, true);
}

void libakaze::bitwise_xor(InputArray a, InputArray b, OutputArray c, InputArray mask)
{
    BinaryFunc f = (BinaryFunc)GET_OPTIMIZED(xor8u);
    binary_op(a, b, c, mask, &f, true);
}



/****************************************************************************************\
*                                      add/subtract                                      *
\****************************************************************************************/

namespace libakaze
{

static int actualScalarDepth(const double* data, int len)
{
    int i = 0, minval = INT_MAX, maxval = INT_MIN;
    for(; i < len; ++i)
    {
        int ival = cvRound(data[i]);
        if( ival != data[i] )
            break;
        minval = MIN(minval, ival);
        maxval = MAX(maxval, ival);
    }
    return i < len ? CV_64F :
        minval >= 0 && maxval <= (int)UCHAR_MAX ? CV_8U :
        minval >= (int)SCHAR_MIN && maxval <= (int)SCHAR_MAX ? CV_8S :
        minval >= 0 && maxval <= (int)USHRT_MAX ? CV_16U :
        minval >= (int)SHRT_MIN && maxval <= (int)SHRT_MAX ? CV_16S :
        CV_32S;
}

static void arithm_op(InputArray _src1, InputArray _src2, OutputArray _dst,
               InputArray _mask, int dtype, BinaryFunc* tab, bool muldiv=false, void* usrdata=0)
{
    int kind1 = _src1.kind(), kind2 = _src2.kind();
    Mat src1 = _src1.getMat(), src2 = _src2.getMat();
    bool haveMask = !_mask.empty();
    bool reallocate = false;

    bool src1Scalar = checkScalar(src1, src2.type(), kind1, kind2);
    bool src2Scalar = checkScalar(src2, src1.type(), kind2, kind1);

    if( (kind1 == kind2 || src1.channels() == 1) && src1.dims <= 2 && src2.dims <= 2 &&
        src1.size() == src2.size() && src1.type() == src2.type() &&
        !haveMask && ((!_dst.fixedType() && (dtype < 0 || CV_MAT_DEPTH(dtype) == src1.depth())) ||
                       (_dst.fixedType() && _dst.type() == _src1.type())) &&
        ((src1Scalar && src2Scalar) || (!src1Scalar && !src2Scalar)) )
    {
        _dst.create(src1.size(), src1.type());
        Mat dst = _dst.getMat();
        Size sz = getContinuousSize(src1, src2, dst, src1.channels());
        tab[src1.depth()](src1.data, src1.step, src2.data, src2.step, dst.data, dst.step, sz, usrdata);
        return;
    }

    bool haveScalar = false, swapped12 = false;
    int depth2 = src2.depth();
    if( src1.size != src2.size || src1.channels() != src2.channels() ||
        (kind1 == _InputArray::MATX && (src1.size() == Size(1,4) || src1.size() == Size(1,1))) ||
        (kind2 == _InputArray::MATX && (src2.size() == Size(1,4) || src2.size() == Size(1,1))) )
    {
        if( checkScalar(src1, src2.type(), kind1, kind2) )
        {
            // src1 is a scalar; swap it with src2
            swap(src1, src2);
            swapped12 = true;
        }
        else if( !checkScalar(src2, src1.type(), kind2, kind1) )
            CV_Error( CV_StsUnmatchedSizes,
                     "The operation is neither 'array op array' (where arrays have the same size and the same number of channels), "
                     "nor 'array op scalar', nor 'scalar op array'" );
        haveScalar = true;
        CV_Assert(src2.type() == CV_64F && (src2.rows == 4 || src2.rows == 1));

        if (!muldiv)
        {
            depth2 = actualScalarDepth(src2.ptr<double>(), src1.channels());
            if( depth2 == CV_64F && (src1.depth() < CV_32S || src1.depth() == CV_32F) )
                depth2 = CV_32F;
        }
        else
            depth2 = CV_64F;
    }

    int cn = src1.channels(), depth1 = src1.depth(), wtype;
    BinaryFunc cvtsrc1 = 0, cvtsrc2 = 0, cvtdst = 0;

    if( dtype < 0 )
    {
        if( _dst.fixedType() )
            dtype = _dst.type();
        else
        {
            if( !haveScalar && src1.type() != src2.type() )
                CV_Error(CV_StsBadArg,
                     "When the input arrays in add/subtract/multiply/divide functions have different types, "
                     "the output array type must be explicitly specified");
            dtype = src1.type();
        }
    }
    dtype = CV_MAT_DEPTH(dtype);

    if( depth1 == depth2 && dtype == depth1 )
        wtype = dtype;
    else if( !muldiv )
    {
        wtype = depth1 <= CV_8S && depth2 <= CV_8S ? CV_16S :
                depth1 <= CV_32S && depth2 <= CV_32S ? CV_32S : std::max(depth1, depth2);
        wtype = std::max(wtype, dtype);

        // when the result of addition should be converted to an integer type,
        // and just one of the input arrays is floating-point, it makes sense to convert that input to integer type before the operation,
        // instead of converting the other input to floating-point and then converting the operation result back to integers.
        if( dtype < CV_32F && (depth1 < CV_32F || depth2 < CV_32F) )
            wtype = CV_32S;
    }
    else
    {
        wtype = std::max(depth1, std::max(depth2, CV_32F));
        wtype = std::max(wtype, dtype);
    }

    cvtsrc1 = depth1 == wtype ? 0 : getConvertFunc(depth1, wtype);
    cvtsrc2 = depth2 == depth1 ? cvtsrc1 : depth2 == wtype ? 0 : getConvertFunc(depth2, wtype);
    cvtdst = dtype == wtype ? 0 : getConvertFunc(wtype, dtype);

    dtype = CV_MAKETYPE(dtype, cn);
    wtype = CV_MAKETYPE(wtype, cn);

    size_t esz1 = src1.elemSize(), esz2 = src2.elemSize();
    size_t dsz = CV_ELEM_SIZE(dtype), wsz = CV_ELEM_SIZE(wtype);
    size_t blocksize0 = (size_t)(BLOCK_SIZE + wsz-1)/wsz;
    BinaryFunc copymask = 0;
    Mat mask;

    if( haveMask )
    {
        mask = _mask.getMat();
        CV_Assert( (mask.type() == CV_8UC1 || mask.type() == CV_8SC1) );
        CV_Assert( mask.size == src1.size );
        copymask = getCopyMaskFunc(dsz);
        Mat tdst = _dst.getMat();
        reallocate = tdst.size != src1.size || tdst.type() != dtype;
    }

    AutoBuffer<uchar> _buf;
    uchar *buf, *maskbuf = 0, *buf1 = 0, *buf2 = 0, *wbuf = 0;
    size_t bufesz = (cvtsrc1 ? wsz : 0) + (cvtsrc2 || haveScalar ? wsz : 0) + (cvtdst ? wsz : 0) + (haveMask ? dsz : 0);

    _dst.create(src1.dims, src1.size, dtype);
    Mat dst = _dst.getMat();

    if( haveMask && reallocate )
        dst = Scalar::all(0);

    BinaryFunc func = tab[CV_MAT_DEPTH(wtype)];

    if( !haveScalar )
    {
        const Mat* arrays[] = { &src1, &src2, &dst, &mask, 0 };
        uchar* ptrs[4];

        NAryMatIterator it(arrays, ptrs);
        size_t total = it.size, blocksize = total;

        if( haveMask || cvtsrc1 || cvtsrc2 || cvtdst )
            blocksize = std::min(blocksize, blocksize0);

        _buf.allocate(bufesz*blocksize + 64);
        buf = _buf;
        if( cvtsrc1 )
            buf1 = buf, buf = alignPtr(buf + blocksize*wsz, 16);
        if( cvtsrc2 )
            buf2 = buf, buf = alignPtr(buf + blocksize*wsz, 16);
        wbuf = maskbuf = buf;
        if( cvtdst )
            buf = alignPtr(buf + blocksize*wsz, 16);
        if( haveMask )
            maskbuf = buf;

        for( size_t i = 0; i < it.nplanes; i++, ++it )
        {
            for( size_t j = 0; j < total; j += blocksize )
            {
                int bsz = (int)MIN(total - j, blocksize);
                Size bszn(bsz*cn, 1);
                const uchar *sptr1 = ptrs[0], *sptr2 = ptrs[1];
                uchar* dptr = ptrs[2];
                if( cvtsrc1 )
                {
                    cvtsrc1( sptr1, 0, 0, 0, buf1, 0, bszn, 0 );
                    sptr1 = buf1;
                }
                if( ptrs[0] == ptrs[1] )
                    sptr2 = sptr1;
                else if( cvtsrc2 )
                {
                    cvtsrc2( sptr2, 0, 0, 0, buf2, 0, bszn, 0 );
                    sptr2 = buf2;
                }

                if( !haveMask && !cvtdst )
                    func( sptr1, 0, sptr2, 0, dptr, 0, bszn, usrdata );
                else
                {
                    func( sptr1, 0, sptr2, 0, wbuf, 0, bszn, usrdata );
                    if( !haveMask )
                        cvtdst( wbuf, 0, 0, 0, dptr, 0, bszn, 0 );
                    else if( !cvtdst )
                    {
                        copymask( wbuf, 0, ptrs[3], 0, dptr, 0, Size(bsz, 1), &dsz );
                        ptrs[3] += bsz;
                    }
                    else
                    {
                        cvtdst( wbuf, 0, 0, 0, maskbuf, 0, bszn, 0 );
                        copymask( maskbuf, 0, ptrs[3], 0, dptr, 0, Size(bsz, 1), &dsz );
                        ptrs[3] += bsz;
                    }
                }
                ptrs[0] += bsz*esz1; ptrs[1] += bsz*esz2; ptrs[2] += bsz*dsz;
            }
        }
    }
    else
    {
        const Mat* arrays[] = { &src1, &dst, &mask, 0 };
        uchar* ptrs[3];

        NAryMatIterator it(arrays, ptrs);
        size_t total = it.size, blocksize = std::min(total, blocksize0);

        _buf.allocate(bufesz*blocksize + 64);
        buf = _buf;
        if( cvtsrc1 )
            buf1 = buf, buf = alignPtr(buf + blocksize*wsz, 16);
        buf2 = buf; buf = alignPtr(buf + blocksize*wsz, 16);
        wbuf = maskbuf = buf;
        if( cvtdst )
            buf = alignPtr(buf + blocksize*wsz, 16);
        if( haveMask )
            maskbuf = buf;

        convertAndUnrollScalar( src2, wtype, buf2, blocksize);

        for( size_t i = 0; i < it.nplanes; i++, ++it )
        {
            for( size_t j = 0; j < total; j += blocksize )
            {
                int bsz = (int)MIN(total - j, blocksize);
                Size bszn(bsz*cn, 1);
                const uchar *sptr1 = ptrs[0];
                const uchar* sptr2 = buf2;
                uchar* dptr = ptrs[1];

                if( cvtsrc1 )
                {
                    cvtsrc1( sptr1, 0, 0, 0, buf1, 0, bszn, 0 );
                    sptr1 = buf1;
                }

                if( swapped12 )
                    std::swap(sptr1, sptr2);

                if( !haveMask && !cvtdst )
                    func( sptr1, 0, sptr2, 0, dptr, 0, bszn, usrdata );
                else
                {
                    func( sptr1, 0, sptr2, 0, wbuf, 0, bszn, usrdata );
                    if( !haveMask )
                        cvtdst( wbuf, 0, 0, 0, dptr, 0, bszn, 0 );
                    else if( !cvtdst )
                    {
                        copymask( wbuf, 0, ptrs[2], 0, dptr, 0, Size(bsz, 1), &dsz );
                        ptrs[2] += bsz;
                    }
                    else
                    {
                        cvtdst( wbuf, 0, 0, 0, maskbuf, 0, bszn, 0 );
                        copymask( maskbuf, 0, ptrs[2], 0, dptr, 0, Size(bsz, 1), &dsz );
                        ptrs[2] += bsz;
                    }
                }
                ptrs[0] += bsz*esz1; ptrs[1] += bsz*dsz;
            }
        }
    }
}

static BinaryFunc* getAddTab()
{
    static BinaryFunc addTab[] =
    {
        (BinaryFunc)GET_OPTIMIZED(add8u), (BinaryFunc)GET_OPTIMIZED(add8s),
        (BinaryFunc)GET_OPTIMIZED(add16u), (BinaryFunc)GET_OPTIMIZED(add16s),
        (BinaryFunc)GET_OPTIMIZED(add32s),
        (BinaryFunc)GET_OPTIMIZED(add32f), (BinaryFunc)add64f,
        0
    };

    return addTab;
}

static BinaryFunc* getSubTab()
{
    static BinaryFunc subTab[] =
    {
        (BinaryFunc)GET_OPTIMIZED(sub8u), (BinaryFunc)GET_OPTIMIZED(sub8s),
        (BinaryFunc)GET_OPTIMIZED(sub16u), (BinaryFunc)GET_OPTIMIZED(sub16s),
        (BinaryFunc)GET_OPTIMIZED(sub32s),
        (BinaryFunc)GET_OPTIMIZED(sub32f), (BinaryFunc)sub64f,
        0
    };

    return subTab;
}



}

void libakaze::add( InputArray src1, InputArray src2, OutputArray dst,
          InputArray mask, int dtype )
{
    arithm_op(src1, src2, dst, mask, dtype, getAddTab() );
}

void libakaze::subtract( InputArray src1, InputArray src2, OutputArray dst,
               InputArray mask, int dtype )
{
#ifdef HAVE_TEGRA_OPTIMIZATION
    if (mask.empty() && src1.depth() == CV_8U && src2.depth() == CV_8U)
    {
        if (dtype == -1 && dst.fixedType())
            dtype = dst.depth();

        if (!dst.fixedType() || dtype == dst.depth())
        {
            if (dtype == CV_16S)
            {
                Mat _dst = dst.getMat();
                if(tegra::subtract_8u8u16s(src1.getMat(), src2.getMat(), _dst))
                    return;
            }
            else if (dtype == CV_32F)
            {
                Mat _dst = dst.getMat();
                if(tegra::subtract_8u8u32f(src1.getMat(), src2.getMat(), _dst))
                    return;
            }
            else if (dtype == CV_8S)
            {
                Mat _dst = dst.getMat();
                if(tegra::subtract_8u8u8s(src1.getMat(), src2.getMat(), _dst))
                    return;
            }
        }
    }
#endif
    arithm_op(src1, src2, dst, mask, dtype, getSubTab() );
}



/****************************************************************************************\
*                                    multiply/divide                                     *
\****************************************************************************************/

namespace libakaze
{

template<typename T, typename WT> static void
mul_( const T* src1, size_t step1, const T* src2, size_t step2,
      T* dst, size_t step, Size size, WT scale )
{
    step1 /= sizeof(src1[0]);
    step2 /= sizeof(src2[0]);
    step /= sizeof(dst[0]);

    if( scale == (WT)1. )
    {
        for( ; size.height--; src1 += step1, src2 += step2, dst += step )
        {
            int i=0;
            #if CV_ENABLE_UNROLLED
            for(; i <= size.width - 4; i += 4 )
            {
                T t0;
                T t1;
                t0 = saturate_cast<T>(src1[i  ] * src2[i  ]);
                t1 = saturate_cast<T>(src1[i+1] * src2[i+1]);
                dst[i  ] = t0;
                dst[i+1] = t1;

                t0 = saturate_cast<T>(src1[i+2] * src2[i+2]);
                t1 = saturate_cast<T>(src1[i+3] * src2[i+3]);
                dst[i+2] = t0;
                dst[i+3] = t1;
            }
            #endif
            for( ; i < size.width; i++ )
                dst[i] = saturate_cast<T>(src1[i] * src2[i]);
        }
    }
    else
    {
        for( ; size.height--; src1 += step1, src2 += step2, dst += step )
        {
            int i = 0;
            #if CV_ENABLE_UNROLLED
            for(; i <= size.width - 4; i += 4 )
            {
                T t0 = saturate_cast<T>(scale*(WT)src1[i]*src2[i]);
                T t1 = saturate_cast<T>(scale*(WT)src1[i+1]*src2[i+1]);
                dst[i] = t0; dst[i+1] = t1;

                t0 = saturate_cast<T>(scale*(WT)src1[i+2]*src2[i+2]);
                t1 = saturate_cast<T>(scale*(WT)src1[i+3]*src2[i+3]);
                dst[i+2] = t0; dst[i+3] = t1;
            }
            #endif
            for( ; i < size.width; i++ )
                dst[i] = saturate_cast<T>(scale*(WT)src1[i]*src2[i]);
        }
    }
}

template<typename T> static void
div_( const T* src1, size_t step1, const T* src2, size_t step2,
      T* dst, size_t step, Size size, double scale )
{
    step1 /= sizeof(src1[0]);
    step2 /= sizeof(src2[0]);
    step /= sizeof(dst[0]);

    for( ; size.height--; src1 += step1, src2 += step2, dst += step )
    {
        int i = 0;
        #if CV_ENABLE_UNROLLED
        for( ; i <= size.width - 4; i += 4 )
        {
            if( src2[i] != 0 && src2[i+1] != 0 && src2[i+2] != 0 && src2[i+3] != 0 )
            {
                double a = (double)src2[i] * src2[i+1];
                double b = (double)src2[i+2] * src2[i+3];
                double d = scale/(a * b);
                b *= d;
                a *= d;

                T z0 = saturate_cast<T>(src2[i+1] * ((double)src1[i] * b));
                T z1 = saturate_cast<T>(src2[i] * ((double)src1[i+1] * b));
                T z2 = saturate_cast<T>(src2[i+3] * ((double)src1[i+2] * a));
                T z3 = saturate_cast<T>(src2[i+2] * ((double)src1[i+3] * a));

                dst[i] = z0; dst[i+1] = z1;
                dst[i+2] = z2; dst[i+3] = z3;
            }
            else
            {
                T z0 = src2[i] != 0 ? saturate_cast<T>(src1[i]*scale/src2[i]) : 0;
                T z1 = src2[i+1] != 0 ? saturate_cast<T>(src1[i+1]*scale/src2[i+1]) : 0;
                T z2 = src2[i+2] != 0 ? saturate_cast<T>(src1[i+2]*scale/src2[i+2]) : 0;
                T z3 = src2[i+3] != 0 ? saturate_cast<T>(src1[i+3]*scale/src2[i+3]) : 0;

                dst[i] = z0; dst[i+1] = z1;
                dst[i+2] = z2; dst[i+3] = z3;
            }
        }
        #endif
        for( ; i < size.width; i++ )
            dst[i] = src2[i] != 0 ? saturate_cast<T>(src1[i]*scale/src2[i]) : 0;
    }
}

template<typename T> static void
recip_( const T*, size_t, const T* src2, size_t step2,
        T* dst, size_t step, Size size, double scale )
{
    step2 /= sizeof(src2[0]);
    step /= sizeof(dst[0]);

    for( ; size.height--; src2 += step2, dst += step )
    {
        int i = 0;
        #if CV_ENABLE_UNROLLED
        for( ; i <= size.width - 4; i += 4 )
        {
            if( src2[i] != 0 && src2[i+1] != 0 && src2[i+2] != 0 && src2[i+3] != 0 )
            {
                double a = (double)src2[i] * src2[i+1];
                double b = (double)src2[i+2] * src2[i+3];
                double d = scale/(a * b);
                b *= d;
                a *= d;

                T z0 = saturate_cast<T>(src2[i+1] * b);
                T z1 = saturate_cast<T>(src2[i] * b);
                T z2 = saturate_cast<T>(src2[i+3] * a);
                T z3 = saturate_cast<T>(src2[i+2] * a);

                dst[i] = z0; dst[i+1] = z1;
                dst[i+2] = z2; dst[i+3] = z3;
            }
            else
            {
                T z0 = src2[i] != 0 ? saturate_cast<T>(scale/src2[i]) : 0;
                T z1 = src2[i+1] != 0 ? saturate_cast<T>(scale/src2[i+1]) : 0;
                T z2 = src2[i+2] != 0 ? saturate_cast<T>(scale/src2[i+2]) : 0;
                T z3 = src2[i+3] != 0 ? saturate_cast<T>(scale/src2[i+3]) : 0;

                dst[i] = z0; dst[i+1] = z1;
                dst[i+2] = z2; dst[i+3] = z3;
            }
        }
        #endif
        for( ; i < size.width; i++ )
            dst[i] = src2[i] != 0 ? saturate_cast<T>(scale/src2[i]) : 0;
    }
}


static void mul8u( const uchar* src1, size_t step1, const uchar* src2, size_t step2,
                   uchar* dst, size_t step, Size sz, void* scale)
{
    mul_(src1, step1, src2, step2, dst, step, sz, (float)*(const double*)scale);
}

static void mul8s( const schar* src1, size_t step1, const schar* src2, size_t step2,
                   schar* dst, size_t step, Size sz, void* scale)
{
    mul_(src1, step1, src2, step2, dst, step, sz, (float)*(const double*)scale);
}

static void mul16u( const ushort* src1, size_t step1, const ushort* src2, size_t step2,
                    ushort* dst, size_t step, Size sz, void* scale)
{
    mul_(src1, step1, src2, step2, dst, step, sz, (float)*(const double*)scale);
}

static void mul16s( const short* src1, size_t step1, const short* src2, size_t step2,
                    short* dst, size_t step, Size sz, void* scale)
{
    mul_(src1, step1, src2, step2, dst, step, sz, (float)*(const double*)scale);
}

static void mul32s( const int* src1, size_t step1, const int* src2, size_t step2,
                    int* dst, size_t step, Size sz, void* scale)
{
    mul_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void mul32f( const float* src1, size_t step1, const float* src2, size_t step2,
                    float* dst, size_t step, Size sz, void* scale)
{
    mul_(src1, step1, src2, step2, dst, step, sz, (float)*(const double*)scale);
}

static void mul64f( const double* src1, size_t step1, const double* src2, size_t step2,
                    double* dst, size_t step, Size sz, void* scale)
{
    mul_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void div8u( const uchar* src1, size_t step1, const uchar* src2, size_t step2,
                   uchar* dst, size_t step, Size sz, void* scale)
{
    if( src1 )
        div_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
    else
        recip_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void div8s( const schar* src1, size_t step1, const schar* src2, size_t step2,
                  schar* dst, size_t step, Size sz, void* scale)
{
    div_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void div16u( const ushort* src1, size_t step1, const ushort* src2, size_t step2,
                    ushort* dst, size_t step, Size sz, void* scale)
{
    div_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void div16s( const short* src1, size_t step1, const short* src2, size_t step2,
                    short* dst, size_t step, Size sz, void* scale)
{
    div_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void div32s( const int* src1, size_t step1, const int* src2, size_t step2,
                    int* dst, size_t step, Size sz, void* scale)
{
    div_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void div32f( const float* src1, size_t step1, const float* src2, size_t step2,
                    float* dst, size_t step, Size sz, void* scale)
{
    div_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void div64f( const double* src1, size_t step1, const double* src2, size_t step2,
                    double* dst, size_t step, Size sz, void* scale)
{
    div_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void recip8u( const uchar* src1, size_t step1, const uchar* src2, size_t step2,
                  uchar* dst, size_t step, Size sz, void* scale)
{
    recip_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void recip8s( const schar* src1, size_t step1, const schar* src2, size_t step2,
                  schar* dst, size_t step, Size sz, void* scale)
{
    recip_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void recip16u( const ushort* src1, size_t step1, const ushort* src2, size_t step2,
                   ushort* dst, size_t step, Size sz, void* scale)
{
    recip_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void recip16s( const short* src1, size_t step1, const short* src2, size_t step2,
                   short* dst, size_t step, Size sz, void* scale)
{
    recip_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void recip32s( const int* src1, size_t step1, const int* src2, size_t step2,
                   int* dst, size_t step, Size sz, void* scale)
{
    recip_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void recip32f( const float* src1, size_t step1, const float* src2, size_t step2,
                   float* dst, size_t step, Size sz, void* scale)
{
    recip_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}

static void recip64f( const double* src1, size_t step1, const double* src2, size_t step2,
                   double* dst, size_t step, Size sz, void* scale)
{
    recip_(src1, step1, src2, step2, dst, step, sz, *(const double*)scale);
}


static BinaryFunc* getMulTab()
{
    static BinaryFunc mulTab[] =
    {
        (BinaryFunc)mul8u, (BinaryFunc)mul8s, (BinaryFunc)mul16u,
        (BinaryFunc)mul16s, (BinaryFunc)mul32s, (BinaryFunc)mul32f,
        (BinaryFunc)mul64f, 0
    };

    return mulTab;
}

static BinaryFunc* getDivTab()
{
    static BinaryFunc divTab[] =
    {
        (BinaryFunc)div8u, (BinaryFunc)div8s, (BinaryFunc)div16u,
        (BinaryFunc)div16s, (BinaryFunc)div32s, (BinaryFunc)div32f,
        (BinaryFunc)div64f, 0
    };

    return divTab;
}

static BinaryFunc* getRecipTab()
{
    static BinaryFunc recipTab[] =
    {
        (BinaryFunc)recip8u, (BinaryFunc)recip8s, (BinaryFunc)recip16u,
        (BinaryFunc)recip16s, (BinaryFunc)recip32s, (BinaryFunc)recip32f,
        (BinaryFunc)recip64f, 0
    };

    return recipTab;
}

}



void libakaze::divide(InputArray src1, InputArray src2,
                OutputArray dst, double scale, int dtype)
{
    arithm_op(src1, src2, dst, noArray(), dtype, getDivTab(), true, &scale);
}

void libakaze::divide(double scale, InputArray src2,
                OutputArray dst, int dtype)
{
    arithm_op(src2, src2, dst, noArray(), dtype, getRecipTab(), true, &scale);
}

/****************************************************************************************\
*                                      addWeighted                                       *
\****************************************************************************************/



/****************************************************************************************\
*                                          compare                                       *
\****************************************************************************************/



/****************************************************************************************\
*                                        inRange                                         *
\****************************************************************************************/


/****************************************************************************************\
*                                Earlier API: cvAdd etc.                                 *
\****************************************************************************************/


/* End of file. */
