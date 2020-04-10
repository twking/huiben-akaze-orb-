#include"cmptkeypoint.h"
#include "fast.h"

#define DBL_EPSILON     2.2204460492503131e-016 /* smallest such that 1.0+DBL_EPSILON != 1.0 */
static int step;
static KEYPOINT_F cmpKeypoint[1000];
static const float atan2_p1 = 0.9997878412794807f*(float)(180/CV_PI);
static const float atan2_p3 = -0.3258083974640975f*(float)(180/CV_PI);
static const float atan2_p5 = 0.1555786518463281f*(float)(180/CV_PI);
static const float atan2_p7 = -0.04432655554792128f*(float)(180/CV_PI);

struct keypointcmp
{
    inline bool operator()(const KEYPOINT_F& kp1, const KEYPOINT_F& kp2) const
    {
		return kp1.val> kp2.val;
    }
};
static int twHarrisResponses(unsigned char *pImag,const int xs,pSKEPOINT pts,int keyPointNums);
static float fastAtan2( float y, float x );
static void IC_Angle(unsigned char *pImag, const int half_k,pSKEPOINT keypoints,const int* u_max);
static void computeOrientation(unsigned char *pImag, pSKEPOINT keypoints,int halfPatchSize, const int* umax,int keynums);


int twcomputeKeyPoints(unsigned char *pImag,const int xs,const int ys,pSKEPOINT pKeyPoint,int maxPointNums,int threshold,int edgeThreshold, int patchSize)
{
	int tempn_points; 
    int halfPatchSize = patchSize / 2;
	int v, v0, vmax = twcvFloor(halfPatchSize * sqrt(2.f) / 2 + 1);
    int vmin = twcvCeil(halfPatchSize * sqrt(2.f) / 2);
    //std::vector<int> umax(halfPatchSize + 2);
	if(halfPatchSize>30)
		return 0;
	int umax[33];
	edgeThreshold += halfPatchSize+1;
	step = xs;
	tempn_points = fast9_16(pImag,xs,ys,cmpKeypoint,2*maxPointNums,threshold,edgeThreshold);
	if(tempn_points==0)
		return 0;
	twHarrisResponses(pImag,xs,cmpKeypoint,tempn_points);
	if(tempn_points>maxPointNums)
	{
		 std::nth_element(cmpKeypoint,cmpKeypoint + maxPointNums, cmpKeypoint+tempn_points, keypointcmp());
		 tempn_points = maxPointNums;
		 memcpy((void*)pKeyPoint,(void *)cmpKeypoint,sizeof(KEYPOINT_F)*maxPointNums);
	}
	else
	{
		memcpy((void*)pKeyPoint,(void *)cmpKeypoint,sizeof(KEYPOINT_F)*tempn_points);
	}
    for (v = 0; v <= vmax; ++v)
        umax[v] = twcvRound(sqrt((float)halfPatchSize * halfPatchSize - v * v));

    // Make sure we are symmetric
    for (v = halfPatchSize, v0 = 0; v >= vmin; --v)
    {
        while (umax[v0] == umax[v0 + 1])
            ++v0;
        umax[v] = v0;
        ++v0;
    }
    computeOrientation(pImag, pKeyPoint, halfPatchSize, umax,tempn_points);
	return tempn_points;
}

static int twHarrisResponses(unsigned char *pImag,const int xs,pSKEPOINT pts,int keyPointNums)
{

	int blockSize = BLOCKSIZE;
	float harris_k = HARRIS_K;
    int ptidx, ptsize = keyPointNums;

    const uchar* ptr00 = pImag;
    int step = xs;
    int r = blockSize/2;

    float scale = (1 << 2) * blockSize * 255.0f;
    scale = 1.0f / scale;
    float scale_sq_sq = scale * scale * scale * scale;

	int ofsbuf[BLOCKSIZE*BLOCKSIZE];// =(int*)twfastMalloc(blockSize*blockSize*sizeof(int),16);
    int* ofs = ofsbuf;
    for( int i = 0; i < blockSize; i++ )
        for( int j = 0; j < blockSize; j++ )
            ofs[i*blockSize + j] = (int)(i*step + j);

    for( ptidx = 0; ptidx < ptsize; ptidx++ )
    {
        int x0 = pts[ptidx].x - r;
        int y0 = pts[ptidx].y - r;

        const uchar* ptr0 = ptr00 + y0*step + x0;
        int a = 0, b = 0, c = 0;

        for( int k = 0; k < blockSize*blockSize; k++ )
        {
            const uchar* ptr = ptr0 + ofs[k];
            int Ix = (ptr[1] - ptr[-1])*2 + (ptr[-step+1] - ptr[-step-1]) + (ptr[step+1] - ptr[step-1]);
            int Iy = (ptr[step] - ptr[-step])*2 + (ptr[step-1] - ptr[-step-1]) + (ptr[step+1] - ptr[-step+1]);
            a += Ix*Ix;
            b += Iy*Iy;
            c += Ix*Iy;
        }
		pts[ptidx].val= ((float)a * b - (float)c * c -harris_k * ((float)a + b) * ((float)a + b))*scale_sq_sq;
    }
	//twfastFree((void *)ofsbuf);
	return 1;
}
static float fastAtan2( float y, float x )
{
    float ax,ay;
	ax= x>0?x:-x;//std::abs(x), 
	ay=y>0?y:-y;// std::abs(y);
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
static void IC_Angle(unsigned char *pImag, const int half_k,pSKEPOINT keypoints,const int*u_max)
{
    int m_01 = 0, m_10 = 0;
    const uchar* center =pImag+step*(keypoints->y)+keypoints->x;

    // Treat the center line differently, v=0
    for (int u = -half_k; u <= half_k; ++u)
        m_10 += u * center[u];

    // Go line by line in the circular patch
    for (int v = 1; v <= half_k; ++v)
    {
        // Proceed over the two lines
        int v_sum = 0;
        int d = u_max[v];
        for (int u = -d; u <= d; ++u)
        {
            int val_plus = center[u + v*step], val_minus = center[u - v*step];
            v_sum += (val_plus - val_minus);
            m_10 += u * (val_plus + val_minus);
        }
        m_01 += v * v_sum;
    }
	keypoints->ang = fastAtan2((float)m_01, (float)m_10);

}
static void computeOrientation(unsigned char *pImag, pSKEPOINT keypoints,int halfPatchSize, const int* umax,int keynums)
{
    // Process each keypoint
    for (int i=0;i<keynums;i++)
    {
		IC_Angle(pImag, halfPatchSize,keypoints+i, umax);
    }
}



