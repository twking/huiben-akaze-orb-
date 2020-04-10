#ifndef _TW_TYPEDEF__
#define _TW_TYPEDEF__

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <iterator>
#include <stdio.h>

using namespace std;
#pragma warning(disable:4996)

#ifdef __cplusplus
extern "C" {
#endif


#define T2X_SYSTEM   'T'
#define GK_SYSTEM	 'G'

#ifndef USING_MODES
#define USING_MODES T2X_SYSTEM
#endif



#define F_MIN(a,b) ((a)>(b)?(b):(a))
#define F_MAX(a,b) ((a)>(b)?(a):(b))


#if USING_MODES == T2X_SYSTEM
#define KEYPOINTNUMS 50 //(50<=KEYPOINTNUMS<=200)
#define DOMEKEYPOINTNUMS 200 //(200<=DOMEKEYPOINTNUMS<=500)
#define VERSION_TYPE "T2X_0"
#elif USING_MODES== GK_SYSTEM
#define KEYPOINTNUMS 100 //(50<=KEYPOINTNUMS<=200)
#define DOMEKEYPOINTNUMS 300 //(200<=DOMEKEYPOINTNUMS<=500)
#define VERSION_TYPE "GK_0"
#endif



typedef unsigned char uchar;
typedef unsigned short ushort;

typedef struct _KEYPOINT
{
	float x;
	float y;
	float val;//�ǵ��ǿ��
}KEYPOINT;
typedef struct _KEYPOINT * pKEYPOINT;

typedef struct _POINYXY
{
	float x;
	float y;
}POINYXY;
typedef struct _POINYXY  Point2f;

typedef struct _DESCRIP 
{
	unsigned char descrip[32];
}DESCRIP;
typedef struct _DESCRIP * pDESCRIP;

typedef struct _TEMPLATEDES 
{
	int pages_key_nums;
	char pages_srial_number[32];//��ҳ�����к�
	void* descrip;//pages descripe DESCRIP
	void* ponit;//pages keypoint POINYXY
}TEMPLATEDES;
typedef struct _TEMPLATEDES * pTEMPLATEDES;

typedef struct _BOOKDESCRIP 
{
	
	 _BOOKDESCRIP()
	 {
		memset(versionType,0,8);
		strcpy(versionType,VERSION_TYPE);
	 }
	 char versionType[8];//��Ӱ汾���ͺ�
	 
	 unsigned int filesize;
	 int nums_desc;//�鱾��ҳ��
	 char srial_number[32];//�鱾�����к�
	 void *pDescrip;//ָ���鱾ģ������ָ��
}BOOKDESCRIP;
typedef struct _BOOKDESCRIP *pBOOKDESCRIP;

typedef struct _Imagetype
{  
    int width;  
    int height;  
    unsigned char* imageData;  
}Imagetype;  
typedef struct _Imagetype *pImagetype;

typedef struct _Matchtype
{  
	int templateIndex;  
    int srcIndex;  
	unsigned int distance;  
}Matchtype;  
typedef struct _Matchtype *pMatchtype;

struct INDENTIFOptions {
	enum
	{
		DATABASE_KEY_MODE_MAX = 1, //�������Ӧֵ�ķ�ʽ��ȡ�ǰ������������
		DATABASE_KEY_MODE_MEAN = 2//�Ծ��ȷֲ��ķ�ʽ��ȡ��������
	};
	enum
	{
		IMAG_TEMPLATE_W = (565 * 2),//һ��ʵ�ʴ���ͼƬ�Ŀ�ȵ�2����Ϊ���ֵ
		DEST_TEMPLATE_H = (346 * 2)//һ��ʵ�ʴ���ͼƬ�ĸ߶ȵ�2����Ϊ���ֵ
	};
	INDENTIFOptions()
		: omax(1)
		, nsublevels(1)
		, img_max_width(IMAG_TEMPLATE_W)
		, img_max_height(DEST_TEMPLATE_H)
		, dthreshold(0.0007f)//[0.005f,0.1f]
		, keymode(DATABASE_KEY_MODE_MEAN)
		, keysize(5)//[5,15],һ������������10*10
		, teamplatekeyPointNums(DOMEKEYPOINTNUMS)//[200,500]��ģ�����ݿ�������������
		, srckeyPointNums(KEYPOINTNUMS)//[50,200]�������ļ������������
		, indentCNNums(8)//[8 15] Ҫ����ƥ�䵽������������
		//, rectificationflag(0)//ͼ���Ƿ������־��0��ʾ��������1��ʾ����
		, dmatchratio(0.8f)//[0.5,0.85] ��ֵ

	{
	}

	int omax;
	int nsublevels;
	int img_max_width;
	int img_max_height;
	float dthreshold;
	int keymode;
	int keysize;
	int teamplatekeyPointNums;
	int srckeyPointNums;
	int indentCNNums;
	//int rectificationflag;
	float dmatchratio;

};

#ifdef __cplusplus
}
#endif

#endif