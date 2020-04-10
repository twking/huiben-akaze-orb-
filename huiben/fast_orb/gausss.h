#ifndef _TW_IIRGAUSSS__
#define _TW_IIRGAUSSS__

#include "tpyedef_orbfast.h"

#ifdef __cplusplus
extern "C" {
#endif

int tw_GaussBlur(unsigned char* src,const int xs,const int ys, float sigma);

#ifdef __cplusplus
}
#endif

#endif

