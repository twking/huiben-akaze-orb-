
#include"fast.h"

static void makeOffsets(int pixel[25], int rowStride)
{
    //�ֱ����������飬���ڱ�ʾpatternSizeΪ16��12��8ʱ��Բ�����ض���Բ�ĵ��������λ��
	int k = 0;
    static const int offsets16[][2] =
    {
        {0,  3}, { 1,  3}, { 2,  2}, { 3,  1}, { 3, 0}, { 3, -1}, { 2, -2}, { 1, -3},
        {0, -3}, {-1, -3}, {-2, -2}, {-3, -1}, {-3, 0}, {-3,  1}, {-2,  2}, {-1,  3}
    };
    //��������ͼ��ÿ�е����ظ������õ�Բ�����صľ�������λ��
    for( ; k < 16; k++ )
        pixel[k] = offsets16[k][0] + offsets16[k][1] * rowStride;
    //����Ҫ�������������أ����Ҫѭ���Ķ��г�һЩֵ
    for( ; k < 25; k++ )
        pixel[k] = pixel[k - 16];
}

static  int cornerScore(const unsigned char* ptr, const int pixel[], int threshold,char condstion)
{
    const int K = 8, N = K*3 + 1;
    //vΪ��ǰ����ֵ
    int k, v = ptr[0];
	int b0,b;
	int a0,a;
    short d[25];
    //���㵱ǰ����ֵ����Բ������ֵ֮��Ĳ�ֵ
    for( k = 0; k < N; k++ )
        d[k] = (short)(v - ptr[pixel[k]]);

    //a0Ϊ��ֵ
	if(condstion==2)
	{
		 a0 = threshold;
		//����ǵ�����2ʱ��������ֵ
		for( k = 0; k < 16; k += 2 )
		{
			//aΪd[k+1]��d[k+2]��d[k+3]�е���Сֵ
			a = F_MIN((int)d[k+1], (int)d[k+2]);
			a = F_MIN(a, (int)d[k+3]);
			//���aС����ֵ���������һ��ѭ��
			if( a <= a0 )
				continue;
			//������ֵ
			//aΪ��d[k+1]��d[k+8]�е���Сֵ
			a = F_MIN(a, (int)d[k+4]);
			a = F_MIN(a, (int)d[k+5]);
			a = F_MIN(a, (int)d[k+6]);
			a = F_MIN(a, (int)d[k+7]);
			a = F_MIN(a, (int)d[k+8]);
			//��d[k]��d[k+9]�е���Сֵ��a0�Ƚϣ��ĸ����ĸ���Ϊ�µ���ֵ
			a0 = F_MAX(a0, F_MIN(a, (int)d[k]));
			a0 = F_MAX(a0, F_MIN(a, (int)d[k+9]));
		}
		threshold = a0;
	}
	if(condstion==1)
	{
		//����ǵ�����1ʱ��������ֵ
		b0 = -threshold;
		for( k = 0; k < 16; k += 2 )
		{
			b = F_MAX((int)d[k+1], (int)d[k+2]);
			b = F_MAX(b, (int)d[k+3]);
			b = F_MAX(b, (int)d[k+4]);
			b = F_MAX(b, (int)d[k+5]);
			if( b >= b0 )
				continue;
			b = F_MAX(b, (int)d[k+6]);
			b = F_MAX(b, (int)d[k+7]);
			b = F_MAX(b, (int)d[k+8]);
 
			b0 = F_MIN(b0, F_MAX(b, (int)d[k]));
			b0 = F_MIN(b0, F_MAX(b, (int)d[k+9]));
		}
		threshold = -b0-1;
	}
    return threshold;
}


void* twfastMalloc(unsigned int size,int alignbytes)
{
	unsigned char** adata;
	unsigned char* udata;
	if(alignbytes%4!=0||alignbytes<4)
		return NULL;
   udata = (unsigned char*)malloc(size + sizeof(void*) + alignbytes);
    if(!udata)
        return NULL;
    adata = (unsigned char**)(((unsigned int)(udata) + alignbytes)&(-alignbytes));
    adata[-1] = udata;
    return adata;
}
 
void twfastFree(void* ptr)
{
    if(ptr)
    {
        unsigned char* udata = ((unsigned char**)ptr)[-1];
        free(udata);
    }
}
int fast9_16(unsigned char *pImag,const int xs,const int ys,pSKEPOINT keypoint,int maxPointNums,int threshold,int edgeThreshold)
{
    //KΪԲ���������صĸ���
    //N����ѭ��Բ�ܵ����ص㣬��ΪҪ��β���ӣ�����NҪ��ʵ��Բ������������K+1��
    const int K = 16/2, N = 16 + K + 1;
	const int cols=xs-2*edgeThreshold;
	const int rows=ys-2*edgeThreshold ;
	pSKEPOINT ptempkey;
	pSKEPOINT pkeya[255];
	int tempthreshold=threshold;
	short* cornerpos;
	short* cpbuf[3];
	unsigned char *ptImag;
	short v,vt,count,x ;
	int keypointNums=0;
	short score,minscore=1255;
	int ncorners;    //��⵽�Ľǵ�����
	int i, j, k, pixel[25];
	// threshold_tabΪ��ֵ�б��ڽ�����ֵ�Ƚϵ�ʱ��ֻ���ñ���
    unsigned char threshold_tab[512];
	unsigned char * _buf;
	unsigned char* buf[3];
	//�õ�buf��ĳ�����飬���ڴ洢��ǰ�еĵ÷ֺ�����ֵV
    unsigned char* curr;
	char d; 
	const unsigned char* tab;
    //�õ�cpbuf��ĳ�����飬���ڴ洢��ǰ�еĽǵ�����λ��

	const unsigned char* prev;
    const unsigned char* pprev;
	const unsigned char* ptr;
	short paary[255];
	memset(paary,0,255*sizeof(short));
	if(maxPointNums>1000||maxPointNums<1)
		return 0;
	if(cols<2*edgeThreshold||rows<2*edgeThreshold)
	{
		return 0;
	}
	else
	{
		ptImag = pImag+edgeThreshold*xs+edgeThreshold;
	}
    makeOffsets(pixel, xs);
	memset(keypoint,0,maxPointNums*sizeof(KEYPOINT_F));
	memset(pkeya,0,255*sizeof(void *));
    /*Ϊ��ֵ�б�ֵ���ñ��Ϊ���Σ���һ�δ�threshold_tab[0]��threshold_tab[255 - threshold]��
	ֵΪ1�����ڸ������ֵ��ʾ����ǵ��ж�����2���ڶ��δ�threshold_tab[255 �C threshold]��threshold_tab[255 + threshold]��
	ֵΪ0�����ڸ������ֵ��ʾ���ǽǵ㣻�����δ�threshold_tab[255 + threshold]��threshold_tab[511]��
	ֵΪ2�����ڸ������ֵ��ʾ����ǵ��ж�����1*/
   // for( i = -255; i <= 255; i++ )
   //     threshold_tab[i+255] = (unsigned char)(i < -tempthreshold ? 1 : i > tempthreshold ? 2 : 0);
	memset(&threshold_tab[0], 1,255 - threshold); 
	memset(&threshold_tab[255 - threshold], 0, threshold*2); 
	memset(&threshold_tab[255 + threshold+1], 2, 255 - threshold+1); 
    //����һ���ڴ�ռ�
    _buf =(unsigned char *)twfastMalloc(((cols+16)*3*(sizeof(short) + sizeof(unsigned char)) + 128),16);
	if(_buf==NULL)
		return 0;
  
    /*buf[0��buf[1]��buf[2]�ֱ��ʾͼ���ǰһ�С���ǰ�кͺ�һ�С���Ϊ�ڷǼ���ֵ���ƵĲ���2�У�
	��Ҫ��3��3�Ľǵ������ڽ��бȽϣ������Ҫ���е�ͼ�����ݡ���Ϊֻ�еõ��˵�ǰ�е����ݣ����Զ�����һ����˵��
	�Ŵչ����������е����ݣ��������ķǼ���ֵ���ƵĽ������һ�����ݵĴ�����*/
    buf[0] = _buf; buf[1] = buf[0] + cols; buf[2] = buf[1] + cols;
    //cpbuf�洢�ǵ������λ�ã�Ҳ����Ҫ�������е�����
    cpbuf[0] = (short*)ALIGNPtr(buf[2] + cols, 4) + 1;
    cpbuf[1] = cpbuf[0] + cols + 1;
    cpbuf[2] = cpbuf[1] + cols + 1;
    memset(buf[0], 0, cols*3);    //buf�����ڴ�����
    //��������ͼ�����أ�Ѱ�ҽǵ�
    //����Բ�İ뾶Ϊ3�����أ����ͼ������ܱ߽綼����3�����صĿ��
    for(i = 3; i < rows-2; i++)
    {
        //�õ�ͼ���е��׵�ַָ��
        ptr =ptImag+i*xs+ 3;
        //�õ�buf��ĳ�����飬���ڴ洢��ǰ�еĵ÷ֺ�����ֵV
        curr = buf[(i - 3)%3];
        //�õ�cpbuf��ĳ�����飬���ڴ洢��ǰ�еĽǵ�����λ��
        cornerpos = cpbuf[(i - 3)%3];
        ncorners = 0;    //��⵽�Ľǵ�����
		memset(curr, 0, cols);    //����
        if( i < rows - 3 )
        {
            //ÿһ�ж�����3�����صĿ��
            j = 3;
            for( ; j < cols - 3; j++, ptr++ )
            {
                //��ǰ���صĻҶ�ֵ
				 v = ptr[0];
                //�ɵ�ǰ���صĻҶ�ֵ��ȷ��������ֵ�б��е�λ��
                tab = &threshold_tab[0] - v + 255;
                //pixel[0]��ʾԲ���ϱ��Ϊ0�����������Բ�������ƫ����
                //ptr[pixel[0]��ʾԲ���ϱ��Ϊ0������ֵ
                //tab[ptr[pixel[0]]]��ʾ����ڵ�ǰ���أ���Բ�ģ�Բ���ϱ��Ϊ0������ֵ����ֵ�б�threshold_tab������ѯ�õ���ֵ�����Ϊ1��˵��I0 < Ip - t�����Ϊ2��˵��I0 > Ip + t�����Ϊ0��˵�� Ip �C t < I0 < Ip + t�����ͨ��tab���Ϳ��Եõ���ǰ�����Ƿ�����ǵ�������
                //���Ϊ0��8����ֱ����Բ���ϵ��������ص㣩���б��е�ֵ����õ�d��d=0˵�����Ϊ0��8��ֵ����0��d=1˵�����Ϊ0��8��ֵ������һ��Ϊ1������һ������Ϊ2��d=2˵�����Ϊ0��8��ֵ������һ��Ϊ2������һ������Ϊ1��d=3˵�����Ϊ0��8��ֵ��һ��Ϊ1����һ��Ϊ2��ֻ�����������������
                d = tab[ptr[pixel[0]]] | tab[ptr[pixel[8]]];
                //d=0˵��Բ���ϲ�����������12����������ǵ���������˵�ǰֵһ�����ǽǵ㣬�����˳��˴�ѭ����������һ��ѭ��
                if( d == 0 )
                    continue;
                //������������ֱ�����������ص���ж�
                d &= tab[ptr[pixel[2]]] | tab[ptr[pixel[10]]];
                d &= tab[ptr[pixel[4]]] | tab[ptr[pixel[12]]];
                d &= tab[ptr[pixel[6]]] | tab[ptr[pixel[14]]];
                //d=0˵������d��������һ��dΪ0�����Կ϶����ǽǵ㣻��һ�������һ��dΪ2������һ��dΪ1�������ҲΪ0����˵��һ��������ǵ�����1������һ������ǵ�����2�����Կ϶�Ҳ����������12����������ͬһ���ǵ������ģ����Ҳһ�����ǽǵ㡣
                if( d == 0 )
                    continue;
                //�����ж�Բ����ʣ������ص�
                d &= tab[ptr[pixel[1]]] | tab[ptr[pixel[9]]];
                d &= tab[ptr[pixel[3]]] | tab[ptr[pixel[11]]];
                d &= tab[ptr[pixel[5]]] | tab[ptr[pixel[13]]];
                d &= tab[ptr[pixel[7]]] | tab[ptr[pixel[15]]];
                //�������if��������˵���п�������ǵ�����2
                if( d & 1 )
                {
                    //vtΪ�����Ľǵ���������Ip �C t��countΪ�������صļ���ֵ
                    vt = v - tempthreshold;
					count = 0;
                    //��������Բ��
                    for( k = 0; k < N; k++ )
                    {
                        x = ptr[pixel[k]];    //��ȡ��Բ���ϵ�����ֵ
                        if(x < vt)    //�����������2
                        {
                            //�������������ж��Ƿ����K��KΪԲ�����ص�һ�룩
                            if( ++count > K )
                            {
                                //�����if��䣬˵���Ѿ��õ�һ���ǵ�
                                //����õ��λ�ã����ѵ�ǰ�еĽǵ�����1
                                cornerpos[ncorners++] = j;
                                 //���зǼ���ֵ���Ƶĵ�һ��������÷ֺ���
								 curr[j] =(unsigned char)cornerScore(ptr,pixel,tempthreshold,2);
                                break;    //�˳�ѭ��
                            }
                        }
                        else
                            count = 0;    //�������صļ���ֵ����
                    }
                }
                //�������if��������˵���п�������ǵ�����1
                if( d & 2 )
                {
                    //vtΪ�����Ľǵ���������Ip + t��countΪ�������صļ���ֵ
                    vt = v + tempthreshold;
					count = 0;
                    //��������Բ��
                    for( k = 0; k < N; k++ )
                    {
                        x = ptr[pixel[k]];    //��ȡ��Բ���ϵ�����ֵ
                        if(x > vt)    //�����������1
                        {
                            //�������������ж��Ƿ����K��KΪԲ�����ص�һ�룩
                            if( ++count > K )
                            {
                                //�����if��䣬˵���Ѿ��õ�һ���ǵ�
                                //����õ��λ�ã����ѵ�ǰ�еĽǵ�����1
                                cornerpos[ncorners++] = j;
                                 //���зǼ���ֵ���Ƶĵ�һ��������÷ֺ���
								 curr[j] =(unsigned char)cornerScore(ptr,pixel,tempthreshold,1);
                                break;    //�˳�ѭ��
                            }
                        }
                        else
                            count = 0;    //�������صļ���ֵ����
                    }
                }
            }
        }
        //���浱ǰ������⵽�Ľǵ���
        cornerpos[-1] = ncorners;
        //i=3˵��ֻ����������һ�е����ݣ������ܽ��зǼ���ֵ���Ƶĵڶ��������Բ������������Ĳ�����ֱ�ӽ�����һ��ѭ��
        if( i == 3 )
            continue;
        //���´����ǽ��зǼ���ֵ���Ƶĵڶ���������3��3�Ľǵ������ڶԵ÷ֺ�����ֵ���зǼ���ֵ���ơ���Ϊ�����������ļ��㣬�Ѿ��õ��˵�ǰ�е����ݣ�
		//���Կ��Խ�����һ�еķǼ���ֵ���ơ��������Ĵ�����е�����һ�еķǼ���ֵ���ơ�
        //��ȡ����һ�к������е�ͼ������
        prev = buf[(i - 4 + 3)%3];
        pprev = buf[(i - 5 + 3)%3];
        //��ȡ����һ������⵽�Ľǵ�λ��
        cornerpos = cpbuf[(i - 4 + 3)%3];
        //��ȡ����һ�еĽǵ���
        ncorners = cornerpos[-1];
        //����һ���ڱ���������⵽�Ľǵ�
        for( k = 0; k < ncorners; k++ )
        {
			
            j = cornerpos[k];    //�õ��ǵ��λ��
            score = prev[j];    //�õ��ýǵ�ĵ÷ֺ���ֵ
            //��3��3�Ľǵ������ڣ����㵱ǰ�ǵ��Ƿ�Ϊ���ֵ���������ѹ������ֵ������
            if( 
				(score > prev[j+1] && score > prev[j-1] &&
                score > pprev[j-1] && score > pprev[j] && score > pprev[j+1] &&
                score > curr[j-1] && score > curr[j] && score > curr[j+1]) )
            {
                //keypoints.push_back(KeyPoint((float)j, (float)(i-1), 7.f, -1, (float)score));
				if(keypointNums==maxPointNums)
				{
					if(score>minscore)
					{
						while(pkeya[minscore]==NULL)
						{
							minscore++;
						}
						paary[score]++;
						paary[minscore]--;
						ptempkey = pkeya[minscore];
						pkeya[minscore] = pkeya[minscore]->pre;
						ptempkey->x = j+edgeThreshold;
						ptempkey->y = i-1+edgeThreshold;
						ptempkey->val = score;
						ptempkey->pre = pkeya[score];
						pkeya[score] = ptempkey;
						tempthreshold = minscore;
					}
					
				}
				else
				{
					keypoint[keypointNums].x=j+edgeThreshold;
					keypoint[keypointNums].y=i-1+edgeThreshold;
					keypoint[keypointNums].val=score;
					keypoint[keypointNums].pre=pkeya[score];
					pkeya[score] = &keypoint[keypointNums];
					keypointNums++;
					if(minscore>score)
						minscore = score;
					paary[score]++;
					if(keypointNums==maxPointNums)
					{
						//tempthreshold = minscore;
						float fk = (float)(rows - i-3.0f)/(i-3.0f);
						int mmmmm=0,iii,ifk;
						ifk=(int)(maxPointNums*fk/(1+fk));
						for(iii=0;iii<255;iii++)
						{
							mmmmm +=paary[iii];
							if(mmmmm>ifk)
								break;
						}
						tempthreshold = (int)(TWK*iii);
						memset(&threshold_tab[255 - tempthreshold], 0, tempthreshold*2); 

					}
					
				}
				
            }
        }
    }
	twfastFree((void *)_buf);
	return keypointNums;
}
