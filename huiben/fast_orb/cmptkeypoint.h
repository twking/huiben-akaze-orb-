#ifndef _TW_CMPTEKEY__
#define _TW_CMPTEKEY__

#include "tpyedef_orbfast.h"

#ifdef __cplusplus
extern "C" {
#endif


int twcomputeKeyPoints(unsigned char *pImag,const int xs,const int ys,pSKEPOINT pKeyPoint,int maxPointNums,int threshold,int edgeThreshold, int patchSize);

#ifdef __cplusplus
}
#endif

#endif