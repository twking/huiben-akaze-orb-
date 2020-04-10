/*
 @function: 1channel 2D gauss filter
 */
#include "gausss.h"

static int twgetGaussianKernel(int *kernel,int n, float sigma)//n=7,sigma=2
{
    const int SMALL_GAUSSIAN_SIZE = 7;
	float cf[14];// = kernel;
    float sigmaX = sigma;// > 0 ? sigma : ((n-1)*0.5 - 1)*0.3 + 0.8;
    float scale2X = -0.5f/(sigmaX*sigmaX);
    float sum = 0;
    int i;
	if(n>14)
		return -1;

    if( n % 2 == 1 && n <= SMALL_GAUSSIAN_SIZE && sigma <= 0)
	{
		return -1;
	}
    for( i = 0; i < n; i++ )
    {
        float x = i - (n-1)*0.5f;
        float t =std::exp(scale2X*x*x);
            cf[i] = (float)t;
            sum += cf[i];
    }
    sum = 1.f/sum;
    for( i = 0; i < n; i++ )
    {
       cf[i] = cf[i]*sum;
	   kernel[i]=(int)(cf[i]*(1<<16));
    }
	return 1;

}


int tw_GaussBlur(unsigned char* src,const int xs,const int ys, float sigma)
{
	int kernel[7];
	int k0,k1,k2,k3,k4,k5,k6;
	int index3,index1,index2;
	int temp=0,n,i;
	unsigned char* pSrc;
	unsigned char* pDst;

	unsigned char*dst= (unsigned char*)malloc(xs*ys*sizeof(char));
    if(dst==NULL)
        return -1;
	if(twgetGaussianKernel(kernel,7,sigma)!=1)
	{
		if(dst!=NULL) free(dst);
		return -1;
	}
	k0=kernel[0];
	k1=kernel[1];	
	k2=kernel[2];
	k3=kernel[3];
	k4=kernel[4];
	k5=kernel[5];
	k6=kernel[6];

	index1=xs;
	index2=2*xs;
	index3=3*xs;

	for(n=3;n<ys-3;n++)
	{
		pSrc = src+n*xs+3;
		pDst = dst+n*xs+3;
		for(i=3;i<xs-3;i++)
		{
			temp=k0**(pSrc-3)+k1**(pSrc-2)+k2**(pSrc-1)+k3**(pSrc)+k4**(pSrc+1)+k5**(pSrc+2)+k6**(pSrc+3)+22768;
			temp = (temp>>16)+1;
			if(temp>255)
			{
				temp=255;
			}
			*pDst  = temp;
			pSrc ++;
			pDst ++;
			
		}
	}
	for(n=3;n<xs-3;n++)
	{
		pSrc = dst+3*xs+n;
		pDst = src+3*xs+n;
		for(i=3;i<ys-3;i++)
		{
			temp=k0**(pSrc-index3)+k1**(pSrc-index2)+k2**(pSrc-index1)+k3**pSrc+k4**(pSrc+index1)+k5**(pSrc+index2)+k6**(pSrc+index3);//+32768;
			temp = (temp>>16)+1;
			if(temp>255)
			{
				temp=255;
			}
			*pDst  = temp;
			pDst +=xs;
			pSrc +=xs;
		}
	}
	if(dst!=NULL) free(dst);
	return 1;

}




