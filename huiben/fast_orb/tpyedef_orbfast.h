#ifndef _ORBFAST_TYPEDEF__
#define _ORBFAST_TYPEDEF__

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <iterator>
#include <stdio.h>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

#define CV_MALLOC_ALIGN 16
#define F_MIN(a,b) ((a)>(b)?(b):(a))
#define F_MAX(a,b) ((a)>(b)?(a):(b))
#define ALIGNPtr( ptr,n)( ( (unsigned int) ptr + n) & -n);
#define TWK 0.8
#define CV_PI   3.1415926535897932384626433832795

#ifndef CV_INLINE
#  if defined __cplusplus
#    define CV_INLINE inline
#  elif defined _MSC_VER
#    define CV_INLINE __inline
#  else
#    define CV_INLINE static
#  endif
#endif /* CV_INLINE */

CV_INLINE  int  twcvRound( double value )
{

    double intpart, fractpart;
    fractpart = modf(value, &intpart);
    if ((fabs(fractpart) != 0.5) || ((((int)intpart) % 2) != 0))
        return (int)(value + (value >= 0 ? 0.5 : -0.5));
    else
        return (int)intpart;

}

CV_INLINE  int  twcvFloor( double value )
{
#if defined _MSC_VER && defined _M_X64 || (defined __GNUC__ && defined __SSE2__ && !defined __APPLE__)
    __m128d t = _mm_set_sd( value );
    int i = _mm_cvtsd_si32(t);
    return i - _mm_movemask_pd(_mm_cmplt_sd(t, _mm_cvtsi32_sd(t,i)));
#elif defined __GNUC__
    int i = (int)value;
    return i - (i > value);
#else
    int i = twcvRound(value);
    float diff = (float)(value - i);
    return i - (diff < 0);
#endif
}

CV_INLINE  int  twcvCeil( double value )
{
#if defined _MSC_VER && defined _M_X64 || (defined __GNUC__ && defined __SSE2__&& !defined __APPLE__)
    __m128d t = _mm_set_sd( value );
    int i = _mm_cvtsd_si32(t);
    return i + _mm_movemask_pd(_mm_cmplt_sd(_mm_cvtsi32_sd(t,i), t));
#elif defined __GNUC__
    int i = (int)value;
    return i + (i < value);
#else
    int i = twcvRound(value);
    float diff = (float)(i - value);
    return i + (diff < 0);
#endif
}

#define HARRIS_K  0.04f
#define DESCRIPTOR_SIZE 32
#define BLOCKSIZE 7

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef struct _SKEYPOINT 
{
	short x;//x坐标
	short y;//y坐标
	float val;//角点的强度
	float ang;//角点的角度
	struct _SKEYPOINT *pre;
}KEYPOINT_F;
typedef struct _SKEYPOINT * pSKEPOINT;

typedef struct _DESCRIPTORS 
{
	unsigned char descrip[32];
	struct _SKEYPOINT *pkeypoint;
}TWDESCRIPTORS;
typedef struct _DESCRIPTORS * pTWDESCRIPTORS;

typedef struct _POINYXYchar
{
	char x;
	char y;
}POINYXY_char;


void*twfastMalloc(unsigned int size,int alignbytes);
void twfastFree(void* ptr);

#ifdef __cplusplus
}
#endif

#endif