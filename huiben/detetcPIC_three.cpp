#include "fast_orb/cmptkeypoint.h"
#include "fast_orb/computeDescriptors.h"
#include "fast_orb/gausss.h"
#include "identifyPictures.h"
#include "dsplib/huiben.h"

#if PICWRITE==1
#include "mkdir_io.h"
#include "c_bmp.h"
#endif


#define DETEC_H 180//120
#define DETEC_KEYPOINT 50
#define RATIO_DETEC 0.6


static DESCRIP srcDesrcrp[DETEC_KEYPOINT];
static int kynums_last=0;
/*
功能：检测当前次图片和上一次图片是否是相同的
参数
1：pImage 指向灰度图片数据
返回值：1：当前次和上一次图片相同
		0：当前次和上一次图片不相同
		-1：函数执行失败
*/
 int detetcPIC_three(Imagetype * pImage)
 {

	 int detec_w,fristpoint_h,fristpoint_w;
	 int cn = 0,i;
	 unsigned char *pdata;
	 unsigned char *temppdata;
	 unsigned char *pTempImag;
	 KEYPOINT_F KeyPoint[DETEC_KEYPOINT];
	 TWDESCRIPTORS pdescriptors[DETEC_KEYPOINT];
	 DESCRIP templateDesrcrp[DETEC_KEYPOINT];
	 Matchtype Ptr[DETEC_KEYPOINT];
	 int kynums=0;

	if(pImage->width<DETEC_H+50) return -1;

	detec_w = DETEC_H*pImage->width/pImage->height;
	fristpoint_h = (int)((pImage->height-DETEC_H)/2);
	fristpoint_w = (int)((pImage->width-detec_w)/2);

	 pdata = (unsigned char *)malloc(detec_w*DETEC_H*sizeof(char));
	 if(pdata == NULL)
	 {
		 return -1;
	 }
	 temppdata = pdata;
	 pTempImag = pImage->imageData + fristpoint_h* pImage->width+fristpoint_w;
	for (int i = 0; i < DETEC_H; i++)
	{		
		memcpy((void*)temppdata, (void *)pTempImag, detec_w);
		temppdata += detec_w;
		pTempImag += pImage->width;

	}

#if PICWRITE==1
	ClImage img_2;
	img_2.channels =1;
	img_2.height=DETEC_H;
	img_2.width=detec_w;
	img_2.imageData=pdata;
	clSaveImage("detetcPIC_three.bmp", &img_2) ;
#endif

	 kynums = twcomputeKeyPoints(pdata,detec_w,DETEC_H,KeyPoint,DETEC_KEYPOINT,20,8, 31);
	if(kynums>10)
	{
		tw_GaussBlur(pdata,detec_w,DETEC_H, 2.0f);
		computeDescriptors(pdata,KeyPoint,pdescriptors,kynums,detec_w);
		for(i=0;i<kynums;i++)
		{
			memcpy(templateDesrcrp[i].descrip,pdescriptors[i].descrip,32);
		}
		cn = hanmmingmatch(templateDesrcrp,kynums,srcDesrcrp,kynums_last,Ptr,RATIO_DETEC);
		kynums_last = kynums;
		memset((void *)srcDesrcrp, 0, sizeof(char) * 32 * DETEC_KEYPOINT);
		memcpy((void*)srcDesrcrp, (void *)templateDesrcrp, sizeof(char) * 32 * kynums);
	}
	if(pdata != NULL)
	{
		free(pdata);
		pdata = NULL;
	}
#if PRINTOUT==1
	printf("detetcPIC_three cn=%d\n",cn);
#endif
	if (cn > kynums / 3)
	{
		return 1;
	}
	else
	{
		return 0;
	}
 }
