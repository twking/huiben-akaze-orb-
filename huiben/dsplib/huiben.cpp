#include "huiben.h"
#include "AKAZEFeatures.h"


#define K_DIS 0.3     //0.1<K_DIS<0.5
#define K_ANGLE 15    //5<K_ANGLE<30   
#define K_SWITCH 0.67 //0.5<K_SWITCH<1
#define CV_PI   3.1415926535897932384626433832795
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

int cmpKeypoint(POINYXY *tempPoint, POINYXY *destPoint, int pointNums)
 {
	 int indexflag;
	 POINYXY temp;
	 float minVal, tempDistance;
	 float x = 0, y = 0;
	 int d_cn, a_cn;
	 int j, i, z;
	 int switcNums=0;
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
	 for (int nnn = 0; nnn < pointNums; nnn++)
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
#if PRINTOUT==1	 	 
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
#endif
	 if (d_cn >= switcNums&&a_cn >= switcNums)
	 {
		 return 1;
	 }
	 else
	 {
		 return -1;
	 }
 }
static bool SortByM3(const libakaze::KeyPoint &v1, const libakaze::KeyPoint &v2)//注意：本函数的参数的类型一定要与vector中元素的类型一致  
{
	return v1.response > v2.response;//降序排列  
}
int detectAndCompute(pImagetype pImage, pKEYPOINT pKey,pDESCRIP pDescrip,unsigned int keynums)
{
	unsigned tempsize;
	libakaze::Mat imageDesc;
	vector<libakaze::KeyPoint>srcKpts;

	if((pImage->imageData==NULL)||(keynums>1000))
	{
		return -1;
	}
	libakaze::Mat tempsrc(pImage->height,pImage->width,CV_8UC1);
	memcpy(tempsrc.data,pImage->imageData,pImage->height*pImage->width*sizeof(char));
	libakaze::AKAZEOptions options;
	options.img_width = tempsrc.cols;
	options.img_height = tempsrc.rows;
	libakaze::AKAZEFeatures impl(options);
	impl.Create_Nonlinear_Scale_Space(tempsrc);
	impl.Feature_Detection(srcKpts);
	if(srcKpts.size()>keynums)
	{
		std::sort(srcKpts.begin(),srcKpts.end(), SortByM3);
		srcKpts.erase(srcKpts.begin()+keynums,srcKpts.end());
	}
	impl.Compute_Descriptors(srcKpts,imageDesc);
	tempsize = srcKpts.size();
	for(unsigned int i=0;i<tempsize;i++)
	{
		pKey->x = srcKpts[i].pt.x;
		pKey->y = srcKpts[i].pt.y;
		pKey->val = srcKpts[i].response;
		pKey++;
	}
	if(srcKpts.size()>0)
	{
		memcpy((void*)pDescrip,(void *)imageDesc.data,32*tempsize);
	}
	return tempsize;

}

#define DETEC_H 180//120
#define DETEC_W (DETEC_H*SCR_W/SCR_H)
#define DETEC_KEYPOINT 50
#define RATIO_DETEC 0.6
#define FRISTPOINT_H (int((360-1.5*DETEC_H)/2))
#define FRISTPOINT_W (int((720-1.5*DETEC_W)/2))
 char detetcPIC_three(unsigned char * pImage)
 {
	 return 0;
 }
 /*
static cv::Mat DetecDesc(DETEC_KEYPOINT, 32, CV_8UC1);
 char detetcPIC_three(unsigned char * pImage)
 {
	 int cn = 0;
	 cv::Ptr<cv::BFMatcher> matcher = cv::BFMatcher::create(cv::NORM_HAMMING);
	 std::vector<vector<DMatch>> mathces;
	 cv::Mat img_1(DETEC_H, DETEC_W, CV_8UC1);
	 //resizePic(pImage, pDesImag);
	 unsigned char *pdata = img_1.data;
	 unsigned char *pTempImag;
	 for (int i = 0; i < DETEC_H; i++)
	 {
		 pTempImag = pImage + 2*(i+ FRISTPOINT_H)* SCR_W+FRISTPOINT_W*2;
		 for (int j = 0; j < DETEC_W; j++)
		 {
			 *pdata++ = pTempImag[j << 1];
		 }
	 }
	 std::vector<KeyPoint> keypoints_dest;
	 cv::Ptr<cv::ORB> orb = cv::ORB::create(DETEC_KEYPOINT, 1.2f, 1, 8, 0, 2, ORB::HARRIS_SCORE, 31);
	 orb->detect(img_1, keypoints_dest);
	 cv::Mat dstDesc;
	 orb->compute(img_1, keypoints_dest, dstDesc);
	 if (keypoints_dest.size()>(DETEC_KEYPOINT/2))
	 {
			matcher->knnMatch(dstDesc, DetecDesc, mathces, 2);
			cn = 0;
			for (int j = 0; j<keypoints_dest.size(); j++)
			{
				if (RATIO_DETEC>mathces[j][0].distance / (mathces[j][1].distance + 0.01))
				{
					cn++;
				}
			}
	 }
	 memset((void *)DetecDesc.data, 0, sizeof(char) * 32 * DETEC_KEYPOINT);
	 memcpy((void*)DetecDesc.data, (void *)dstDesc.data, sizeof(char) * 32 * keypoints_dest.size());
	 printf("detetcPIC_three cn=%d\n",cn);
	 if (cn > keypoints_dest.size() / 3)
	 {
		 return 1;
	 }
	 else
	 {
		 return 0;
	 }
 }
 */
