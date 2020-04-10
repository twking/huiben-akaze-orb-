#include "pc_huiben.h"
#pragma warning(disable:4996)

static INDENTIFOptions identifyoptions;

#define K_DIS 0.3     //0.1<K_DIS<0.5
#define K_ANGLE 15    //5<K_ANGLE<30   
#define K_SWITCH 0.67 //0.5<K_SWITCH<1
#define CV_PI   3.1415926535897932384626433832795
#define DBL_EPSILON     2.2204460492503131e-016 /* smallest such that 1.0+DBL_EPSILON != 1.0 */

static char retStr[32 + 1] = { 0 };//sizeof(BOOKDESCRIP::srial_number)
static bool SortByM4(const Matchtype &v1, const Matchtype &v2)//注意：本函数的参数的类型一定要与vector中元素的类型一致  
{
	return v1.distance < v2.distance;//升序排列  
}
static bool SortByM1(const float &v1, const float &v2)//注意：本函数的参数的类型一定要与vector中元素的类型一致  
{
	 return v1 < v2;//升序排列  
} 
static const float atan2_p1 = 0.9997878412794807f*(float)(180 / CV_PI);
static const float atan2_p3 = -0.3258083974640975f*(float)(180 / CV_PI);
static const float atan2_p5 = 0.1555786518463281f*(float)(180 / CV_PI);
static const float atan2_p7 = -0.04432655554792128f*(float)(180 / CV_PI);
static inline float twfastAtan2( float y, float x )
{
    float ax = std::abs(x), ay = std::abs(y);
    float a, c, c2;
    if( ax >= ay )
    {
        c = ay/(ax + (float)DBL_EPSILON);
        c2 = c*c;
        a = (((atan2_p7*c2 + atan2_p5)*c2 + atan2_p3)*c2 + atan2_p1)*c;
    }
    else
    {
        c = ax/(ay + (float)DBL_EPSILON);
        c2 = c*c;
        a = 90.f - (((atan2_p7*c2 + atan2_p5)*c2 + atan2_p3)*c2 + atan2_p1)*c;
    }
    if( x < 0 )
        a = 180.f - a;
    if( y < 0 )
        a = 360.f - a;
    return a;
}
 int cmpKeypoint(POINYXY *tempPoint, POINYXY *destPoint, unsigned int pointNums)
 {
	 unsigned int indexflag;
	 POINYXY temp;
	 float minVal, tempDistance;
	 float x = 0, y = 0;
	 unsigned int d_cn, a_cn;
	 unsigned int j, i, z;
	 unsigned int switcNums=0;
	 vector<float> tempPointDistance, destPointDsistance, k_Distance;
	 vector<float> tempPointAngle, destPointAngle, abs_angle;
	 vector<int> index;
	 if (pointNums<6)
		 return -1;
	 index.push_back(0);
	 for (j = 1; j<pointNums; j++)
	 {

		 x = tempPoint[j - 1].x;
		 y = tempPoint[j - 1].y;
		 minVal = 88888888888888888888.8f;
		 indexflag = j;
		 for (i = j; i<pointNums; i++)
		 {
			 tempDistance = (tempPoint[i].x - x)*(tempPoint[i].x - x) + (tempPoint[i].y - y)*(tempPoint[i].y - y);
			 if (tempDistance<minVal)
			 {
				 minVal = tempDistance;
				 indexflag = i;
			 }
		 }
		 if (indexflag != j)
		 {
			 temp = tempPoint[j];
			 tempPoint[j] = tempPoint[indexflag];
			 tempPoint[indexflag] = temp;

			 temp = destPoint[j];
			 destPoint[j] = destPoint[indexflag];
			 destPoint[indexflag] = temp;
		 }
		 index.push_back(j);
	 }
	 indexflag = pointNums >> 1;
	 for (j = 0, i = 0; i<pointNums - 2; i += 2)
	 {
		 index[i] = indexflag + j;
		 index[i + 1] = j;
		 j++;
	 }
	 int ccccn = 0;
	 for (unsigned int nnn = 0; nnn < pointNums; nnn++)
	 {
		 for (j = nnn + 3; j<pointNums; j++)
		 {

			 i = nnn;
			 if (j == nnn)
			 {
				 continue;
			 }
			 z = j;
			 tempPointDistance.push_back(sqrt((tempPoint[i].x - tempPoint[z].x)*(tempPoint[i].x - tempPoint[z].x) + (tempPoint[i].y - tempPoint[z].y)*(tempPoint[i].y - tempPoint[z].y)));
			 destPointDsistance.push_back(sqrt((destPoint[i].x - destPoint[z].x)*(destPoint[i].x - destPoint[z].x) + (destPoint[i].y - destPoint[z].y)*(destPoint[i].y - destPoint[z].y)));
			 tempPointAngle.push_back(twfastAtan2(tempPoint[i].x - tempPoint[z].x, tempPoint[i].y - tempPoint[z].y));
			 destPointAngle.push_back(twfastAtan2(destPoint[i].x - destPoint[z].x, destPoint[i].y - destPoint[z].y));
			 abs_angle.push_back(abs(tempPointAngle[ccccn] - destPointAngle[ccccn]));
			 if (abs_angle[ccccn]>180)
				 abs_angle[ccccn] = 360 - abs_angle[ccccn];
			 destPointAngle[ccccn] = tempPointAngle[ccccn] + abs_angle[ccccn];
			 k_Distance.push_back(F_MIN(F_MAX(destPointDsistance[ccccn] / (tempPointDistance[ccccn] + 1.0f), 0.01f), 10));
			 ccccn++;
		 }
	 }

	 std::sort(k_Distance.begin(), k_Distance.end(), SortByM1);
	 std::sort(abs_angle.begin(), abs_angle.end(), SortByM1);
	 switcNums = (int)((ccccn)*K_SWITCH + 0.5);
	 d_cn = 1;
	 for (j = 0; j<k_Distance.size(); j++)
	 {
		 z = 1;
		 for (i = j + 1; i<k_Distance.size(); i++)
		 {
			 if ((k_Distance[i] - k_Distance[j])<K_DIS*k_Distance[j] && k_Distance[i]> 0.16&&k_Distance[i]<6)
			 {
				 z++;
				 if (d_cn<z)
					 d_cn = z;
			 }
			 else
			 {
				 break;
			 }
		 }
	 }
	 a_cn = 0;
	 for (j = 0; j<abs_angle.size(); j++)
	 {
		 z = 1;
		 for (i = j + 1; i<abs_angle.size(); i++)
		 {
			 if ((abs_angle[i] - abs_angle[j])<K_ANGLE)
			 {
				 z++;
				 if (a_cn<z)
					 a_cn = z;
			 }
			 else
			 {
				 break;
			 }
		 }
	 }
	 
	 printf("***********switcNums=%d*****************\n", switcNums);
	 printf("k_Distance:d_cn=%d\n", d_cn);
	 for (j = 0; j<k_Distance.size(); j++)
	 {
		 printf("[%d]=%3.2f,", j, k_Distance[j]);
		 if (j % 5 == 0 && j>0)
		 {
			 printf("\n");
		 }
	 }
	 printf("\nabs_angle:a_cn=%d\n", a_cn);
	 for (j = 0; j<abs_angle.size(); j++)
	 {
		 printf("[%d]=%3.1f,", j, abs_angle[j]);
		 if (j % 5 == 0 && j>0)
		 {
			 printf("\n");
		 }
	 }
	 printf("\n");
	 
	 if (d_cn >= switcNums&&a_cn >= switcNums)
	 {
		 return 1;
	 }
	 else
	 {
		 return -1;
	 }
 }
char *pc_identifycontest(void * pdecvToPCDesribe, pBOOKDESCRIP pCONTENTESDescribe)
{
	char * returnNums = NULL;
	DESCRIP *tempsDesc = NULL;
	POINYXY *tempskey = NULL;
	int cn = 0, i, j,retcn;
	int pagesKeyNums;
	unsigned char *preturnData;
	unsigned int addjiaoyan=0,addjiaoyan1=0;
	if (pdecvToPCDesribe == NULL || pCONTENTESDescribe == NULL)
		return NULL;

	pTEMPLATEDES pBOOKDescribe= (pTEMPLATEDES)pdecvToPCDesribe;

	if(pBOOKDescribe->pages_key_nums<10||pBOOKDescribe->pages_key_nums>200) return NULL;

	retcn = sizeof(int) + 32 * sizeof(char) + pBOOKDescribe->pages_key_nums * sizeof(DESCRIP) + pBOOKDescribe->pages_key_nums * sizeof(POINYXY)+sizeof(addjiaoyan);
	preturnData= (unsigned char *)pdecvToPCDesribe;

	if(strcmp(pCONTENTESDescribe->versionType,pBOOKDescribe->pages_srial_number)!=0)
	{
		printf("error:pc_identifycontest:Library version mismatch:%s!=%s\n",pCONTENTESDescribe->versionType,pBOOKDescribe->pages_srial_number);
		return NULL;
	}
	for(i=0;i<retcn-4;i++)
	{
		addjiaoyan +=preturnData[i];
	}
	addjiaoyan1 = preturnData[i]+preturnData[i+1]*(1<<8)+preturnData[i+2]*(1<<16)+preturnData[i+3]*(1<<24);
	if(addjiaoyan1!=addjiaoyan) return NULL;

	DESCRIP *dstDesc = (DESCRIP *)((char *)pdecvToPCDesribe + sizeof(pBOOKDescribe->pages_key_nums) + sizeof(pBOOKDescribe->pages_srial_number));
	POINYXY *keypoints_dest = (POINYXY *)((char *)pdecvToPCDesribe + sizeof(pBOOKDescribe->pages_key_nums) + sizeof(pBOOKDescribe->pages_srial_number) + pBOOKDescribe->pages_key_nums * sizeof(DESCRIP));

	pTEMPLATEDES ptempDesc = (pTEMPLATEDES)pCONTENTESDescribe->pDescrip;
	int pages = pCONTENTESDescribe->nums_desc;

	POINYXY *pScr = new POINYXY[pBOOKDescribe->pages_key_nums];
	POINYXY *pDest = new POINYXY[pBOOKDescribe->pages_key_nums];
	Matchtype *matchval = new Matchtype[pBOOKDescribe->pages_key_nums];
	DESCRIP **IndixtempsDesc = new DESCRIP*[pages];
	POINYXY **IndixtempsKeyc = new POINYXY*[pages];
	int *pkeyNums = new int[pages];

	if (pScr == NULL || pDest == NULL ||matchval == NULL
		||IndixtempsDesc == NULL || IndixtempsKeyc == NULL 
		||pkeyNums == NULL) goto retn;

	if (pCONTENTESDescribe->pDescrip != NULL)
	{
		char *p = (char *)ptempDesc;
		for (i = 0; i<pages; i++)
		{
			pagesKeyNums = ptempDesc->pages_key_nums;
			pkeyNums[i] = pagesKeyNums;
			p += sizeof(ptempDesc->pages_key_nums) + sizeof(ptempDesc->pages_srial_number);//36
			IndixtempsDesc[i] = (DESCRIP *)p;
			p += pagesKeyNums * sizeof(DESCRIP);
			IndixtempsKeyc[i] = (POINYXY *)p;
			p += pagesKeyNums * sizeof(POINYXY);
			ptempDesc = (pTEMPLATEDES)p;
		}
		for (i = 0; i<pages; i++)
		{
			cn = hanmmingmatch(IndixtempsDesc[i], pkeyNums[i], dstDesc, pBOOKDescribe->pages_key_nums, matchval, identifyoptions.dmatchratio);
			printf("contents:In pages %d find Similar characteristics is %d;\n", i, cn);
			if (cn >= identifyoptions.indentCNNums)
			{
				std::sort(matchval, matchval + cn, SortByM4);
				for (j = 0; j<identifyoptions.indentCNNums; j++)
				{
					pDest[j].x = keypoints_dest[matchval[j].srcIndex].x;
					pDest[j].y = keypoints_dest[matchval[j].srcIndex].y;
					pScr[j].x = IndixtempsKeyc[i][matchval[j].templateIndex].x;
					pScr[j].y = IndixtempsKeyc[i][matchval[j].templateIndex].y;

				}
				if (cmpKeypoint(pDest, pScr, identifyoptions.indentCNNums) == 1)
				{
					memset(retStr, 0, sizeof(pCONTENTESDescribe->srial_number) + 1);
					memcpy(retStr, ((char *)IndixtempsDesc[i]) - sizeof(pCONTENTESDescribe->srial_number), sizeof(pCONTENTESDescribe->srial_number));
					returnNums = retStr;
					goto retn;
				}
			}
		}
	}
retn:
	if (pScr != NULL) delete[] pScr;
	if (pDest != NULL) delete[] pDest;
	if (matchval != NULL)delete[] matchval;
	if (IndixtempsDesc != NULL) delete[] IndixtempsDesc;
	if (IndixtempsKeyc != NULL) delete[] IndixtempsKeyc;
	if (pkeyNums != NULL) delete[] pkeyNums;

	return returnNums;
}
int readDesripFromfile(string str, pBOOKDESCRIP pBooksDesc)
{
	FILE *fid;
	fid = fopen(str.c_str(), "rb");
	if (fid == NULL)
	{
		printf("read file failed!\n");
		return 0;
	}
	//获取文件大小
	long tempSize;
	fseek(fid, 0, SEEK_END);
	long lSize = ftell(fid);
	rewind(fid);
	fread(&pBooksDesc->versionType, sizeof(pBooksDesc->versionType)+sizeof(pBooksDesc->filesize)+sizeof(pBooksDesc->nums_desc)+sizeof(pBooksDesc->srial_number), 1, fid);
	tempSize = pBooksDesc->filesize;
	if(strcmp(pBooksDesc->versionType,VERSION_TYPE)!=0)
	{
		printf("Library version mismatch:%s\n",pBooksDesc->versionType);
		return 0;
	}
	if (lSize != tempSize + sizeof(pBooksDesc->versionType)+sizeof(pBooksDesc->filesize)+sizeof(pBooksDesc->nums_desc)+sizeof(pBooksDesc->srial_number))
	{
		printf("file is damage off!\n");
		return 0;
	}
	printf("versionType=%s,lSize=%d,tempSize=%d\n", pBooksDesc->versionType,lSize, tempSize);
	//开辟存储空间
	pBooksDesc->pDescrip = malloc(tempSize);
	if (pBooksDesc->pDescrip == NULL)
	{
		printf("malloc memory spach failed!\n");

		return 0;
	}
	fread(pBooksDesc->pDescrip, sizeof(char), tempSize, fid);
	fclose(fid);
	return lSize;
}

