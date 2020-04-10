#ifndef _TW_CMPTEDSP__
#define _TW_CMPTEDSP__
#include "tpyedef_orbfast.h"
#ifdef __cplusplus
extern "C" {
#endif

void computeDescriptors(const unsigned char *image,pSKEPOINT keypoints, pTWDESCRIPTORS pdescriptors, int nkeyPoint,int step);
	
#ifdef __cplusplus
}
#endif

#endif