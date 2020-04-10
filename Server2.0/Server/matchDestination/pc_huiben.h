#ifndef _TW_HUIBEN_H__
#define _TW_HUIBEN_H__

#include "../pc_dsplib/tw_tpyedef.h"

#ifdef __cplusplus
extern "C" {
#endif

unsigned int hanmmingmatch(pDESCRIP templateDesrcrp,unsigned int templateDesrcrpnums,pDESCRIP srcDesrcrp,unsigned int srcnums,pMatchtype Ptr,float k);
 int cmpKeypoint(POINYXY *tempPoint, POINYXY *destPoint, unsigned int pointNums);
char *pc_identifycontest(void * pdecvToPCDesribe, pBOOKDESCRIP pCONTENTESDescribe);
int readDesripFromfile(string str, pBOOKDESCRIP pBooksDesc);
#ifdef __cplusplus
}
#endif

#endif
