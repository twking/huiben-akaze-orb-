#ifndef _TW_FAST__
#define _TW_FAST__

#include "tpyedef_orbfast.h"

#ifdef __cplusplus
extern "C" {
#endif

int fast9_16(unsigned char *pImag,const int xs,const int ys,pSKEPOINT pKeyPoint,int maxPointNums,int threshold,int edgeThreshold);


#ifdef __cplusplus
}
#endif

#endif