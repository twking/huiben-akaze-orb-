#ifndef _TW_HUIBEN_H__
#define _TW_HUIBEN_H__

#include "tw_tpyedef.h"

#ifdef __cplusplus
extern "C" {
#endif

unsigned int hanmmingmatch(pDESCRIP templateDesrcrp,unsigned int templateDesrcrpnums,pDESCRIP srcDesrcrp,unsigned int srcnums,pMatchtype Ptr,float k);
int detectAndCompute(pImagetype pImage, pKEYPOINT pKey,pDESCRIP pDescrip,unsigned int keynums);
int cmpKeypoint(POINYXY *tempPoint, POINYXY *destPoint, int pointNums); 

#ifdef __cplusplus
}
#endif

#endif
