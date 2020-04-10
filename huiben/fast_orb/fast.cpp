
#include"fast.h"

static void makeOffsets(int pixel[25], int rowStride)
{
    //分别定义三个数组，用于表示patternSize为16，12和8时，圆周像素对于圆心的相对坐标位置
	int k = 0;
    static const int offsets16[][2] =
    {
        {0,  3}, { 1,  3}, { 2,  2}, { 3,  1}, { 3, 0}, { 3, -1}, { 2, -2}, { 1, -3},
        {0, -3}, {-1, -3}, {-2, -2}, {-3, -1}, {-3, 0}, {-3,  1}, {-2,  2}, {-1,  3}
    };
    //代入输入图像每行的像素个数，得到圆周像素的绝对坐标位置
    for( ; k < 16; k++ )
        pixel[k] = offsets16[k][0] + offsets16[k][1] * rowStride;
    //由于要计算连续的像素，因此要循环的多列出一些值
    for( ; k < 25; k++ )
        pixel[k] = pixel[k - 16];
}

static  int cornerScore(const unsigned char* ptr, const int pixel[], int threshold,char condstion)
{
    const int K = 8, N = K*3 + 1;
    //v为当前像素值
    int k, v = ptr[0];
	int b0,b;
	int a0,a;
    short d[25];
    //计算当前像素值与其圆周像素值之间的差值
    for( k = 0; k < N; k++ )
        d[k] = (short)(v - ptr[pixel[k]]);

    //a0为阈值
	if(condstion==2)
	{
		 a0 = threshold;
		//满足角点条件2时，更新阈值
		for( k = 0; k < 16; k += 2 )
		{
			//a为d[k+1]，d[k+2]和d[k+3]中的最小值
			a = F_MIN((int)d[k+1], (int)d[k+2]);
			a = F_MIN(a, (int)d[k+3]);
			//如果a小于阈值，则进行下一次循环
			if( a <= a0 )
				continue;
			//更新阈值
			//a为从d[k+1]到d[k+8]中的最小值
			a = F_MIN(a, (int)d[k+4]);
			a = F_MIN(a, (int)d[k+5]);
			a = F_MIN(a, (int)d[k+6]);
			a = F_MIN(a, (int)d[k+7]);
			a = F_MIN(a, (int)d[k+8]);
			//从d[k]到d[k+9]中的最小值与a0比较，哪个大，哪个作为新的阈值
			a0 = F_MAX(a0, F_MIN(a, (int)d[k]));
			a0 = F_MAX(a0, F_MIN(a, (int)d[k+9]));
		}
		threshold = a0;
	}
	if(condstion==1)
	{
		//满足角点条件1时，更新阈值
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
    //K为圆周连续像素的个数
    //N用于循环圆周的像素点，因为要首尾连接，所以N要比实际圆周像素数量多K+1个
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
	int ncorners;    //检测到的角点数量
	int i, j, k, pixel[25];
	// threshold_tab为阈值列表，在进行阈值比较的时候，只需查该表即可
    unsigned char threshold_tab[512];
	unsigned char * _buf;
	unsigned char* buf[3];
	//得到buf的某个数组，用于存储当前行的得分函数的值V
    unsigned char* curr;
	char d; 
	const unsigned char* tab;
    //得到cpbuf的某个数组，用于存储当前行的角点坐标位置

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
    /*为阈值列表赋值，该表分为三段：第一段从threshold_tab[0]至threshold_tab[255 - threshold]，
	值为1，落在该区域的值表示满足角点判断条件2；第二段从threshold_tab[255 – threshold]至threshold_tab[255 + threshold]，
	值为0，落在该区域的值表示不是角点；第三段从threshold_tab[255 + threshold]至threshold_tab[511]，
	值为2，落在该区域的值表示满足角点判断条件1*/
   // for( i = -255; i <= 255; i++ )
   //     threshold_tab[i+255] = (unsigned char)(i < -tempthreshold ? 1 : i > tempthreshold ? 2 : 0);
	memset(&threshold_tab[0], 1,255 - threshold); 
	memset(&threshold_tab[255 - threshold], 0, threshold*2); 
	memset(&threshold_tab[255 + threshold+1], 2, 255 - threshold+1); 
    //开辟一段内存空间
    _buf =(unsigned char *)twfastMalloc(((cols+16)*3*(sizeof(short) + sizeof(unsigned char)) + 128),16);
	if(_buf==NULL)
		return 0;
  
    /*buf[0、buf[1]和buf[2]分别表示图像的前一行、当前行和后一行。因为在非极大值抑制的步骤2中，
	是要在3×3的角点邻域内进行比较，因此需要三行的图像数据。因为只有得到了当前行的数据，所以对于上一行来说，
	才凑够了连续三行的数据，因此输出的非极大值抑制的结果是上一行数据的处理结果*/
    buf[0] = _buf; buf[1] = buf[0] + cols; buf[2] = buf[1] + cols;
    //cpbuf存储角点的坐标位置，也是需要连续三行的数据
    cpbuf[0] = (short*)ALIGNPtr(buf[2] + cols, 4) + 1;
    cpbuf[1] = cpbuf[0] + cols + 1;
    cpbuf[2] = cpbuf[1] + cols + 1;
    memset(buf[0], 0, cols*3);    //buf数组内存清零
    //遍历整幅图像像素，寻找角点
    //由于圆的半径为3个像素，因此图像的四周边界都留出3个像素的宽度
    for(i = 3; i < rows-2; i++)
    {
        //得到图像行的首地址指针
        ptr =ptImag+i*xs+ 3;
        //得到buf的某个数组，用于存储当前行的得分函数的值V
        curr = buf[(i - 3)%3];
        //得到cpbuf的某个数组，用于存储当前行的角点坐标位置
        cornerpos = cpbuf[(i - 3)%3];
        ncorners = 0;    //检测到的角点数量
		memset(curr, 0, cols);    //清零
        if( i < rows - 3 )
        {
            //每一行都留出3个像素的宽度
            j = 3;
            for( ; j < cols - 3; j++, ptr++ )
            {
                //当前像素的灰度值
				 v = ptr[0];
                //由当前像素的灰度值，确定其在阈值列表中的位置
                tab = &threshold_tab[0] - v + 255;
                //pixel[0]表示圆周上编号为0的像素相对于圆心坐标的偏移量
                //ptr[pixel[0]表示圆周上编号为0的像素值
                //tab[ptr[pixel[0]]]表示相对于当前像素（即圆心）圆周上编号为0的像素值在阈值列表threshold_tab中所查询得到的值，如果为1，说明I0 < Ip - t，如果为2，说明I0 > Ip + t，如果为0，说明 Ip – t < I0 < Ip + t。因此通过tab，就可以得到当前像素是否满足角点条件。
                //编号为0和8（即直径在圆周上的两个像素点）在列表中的值相或后得到d。d=0说明编号为0和8的值都是0；d=1说明编号为0和8的值至少有一个为1，而另一个不能为2；d=2说明编号为0和8的值至少有一个为2，而另一个不能为1；d=3说明编号为0和8的值有一个为1，另一个为2。只可能有这四种情况。
                d = tab[ptr[pixel[0]]] | tab[ptr[pixel[8]]];
                //d=0说明圆周上不可能有连续12个像素满足角点条件，因此当前值一定不是角点，所以退出此次循环，进入下一次循环
                if( d == 0 )
                    continue;
                //继续进行其他直径上两个像素点的判断
                d &= tab[ptr[pixel[2]]] | tab[ptr[pixel[10]]];
                d &= tab[ptr[pixel[4]]] | tab[ptr[pixel[12]]];
                d &= tab[ptr[pixel[6]]] | tab[ptr[pixel[14]]];
                //d=0说明上述d中至少有一个d为0，所以肯定不是角点；另一种情况是一个d为2，而另一个d为1，相与后也为0，这说明一个是满足角点条件1，而另一个满足角点条件2，所以肯定也不会有连续12个像素满足同一个角点条件的，因此也一定不是角点。
                if( d == 0 )
                    continue;
                //继续判断圆周上剩余的像素点
                d &= tab[ptr[pixel[1]]] | tab[ptr[pixel[9]]];
                d &= tab[ptr[pixel[3]]] | tab[ptr[pixel[11]]];
                d &= tab[ptr[pixel[5]]] | tab[ptr[pixel[13]]];
                d &= tab[ptr[pixel[7]]] | tab[ptr[pixel[15]]];
                //如果满足if条件，则说明有可能满足角点条件2
                if( d & 1 )
                {
                    //vt为真正的角点条件，即Ip – t，count为连续像素的计数值
                    vt = v - tempthreshold;
					count = 0;
                    //遍历整个圆周
                    for( k = 0; k < N; k++ )
                    {
                        x = ptr[pixel[k]];    //提取出圆周上的像素值
                        if(x < vt)    //如果满足条件2
                        {
                            //连续计数，并判断是否大于K（K为圆周像素的一半）
                            if( ++count > K )
                            {
                                //进入该if语句，说明已经得到一个角点
                                //保存该点的位置，并把当前行的角点数加1
                                cornerpos[ncorners++] = j;
                                 //进行非极大值抑制的第一步，计算得分函数
								 curr[j] =(unsigned char)cornerScore(ptr,pixel,tempthreshold,2);
                                break;    //退出循环
                            }
                        }
                        else
                            count = 0;    //连续像素的计数值清零
                    }
                }
                //如果满足if条件，则说明有可能满足角点条件1
                if( d & 2 )
                {
                    //vt为真正的角点条件，即Ip + t，count为连续像素的计数值
                    vt = v + tempthreshold;
					count = 0;
                    //遍历整个圆周
                    for( k = 0; k < N; k++ )
                    {
                        x = ptr[pixel[k]];    //提取出圆周上的像素值
                        if(x > vt)    //如果满足条件1
                        {
                            //连续计数，并判断是否大于K（K为圆周像素的一半）
                            if( ++count > K )
                            {
                                //进入该if语句，说明已经得到一个角点
                                //保存该点的位置，并把当前行的角点数加1
                                cornerpos[ncorners++] = j;
                                 //进行非极大值抑制的第一步，计算得分函数
								 curr[j] =(unsigned char)cornerScore(ptr,pixel,tempthreshold,1);
                                break;    //退出循环
                            }
                        }
                        else
                            count = 0;    //连续像素的计数值清零
                    }
                }
            }
        }
        //保存当前行所检测到的角点数
        cornerpos[-1] = ncorners;
        //i=3说明只仅仅计算了一行的数据，还不能进行非极大值抑制的第二步，所以不进行下面代码的操作，直接进入下一次循环
        if( i == 3 )
            continue;
        //以下代码是进行非极大值抑制的第二步，即在3×3的角点邻域内对得分函数的值进行非极大值抑制。因为经过上面代码的计算，已经得到了当前行的数据，
		//所以可以进行上一行的非极大值抑制。因此下面的代码进行的是上一行的非极大值抑制。
        //提取出上一行和上两行的图像像素
        prev = buf[(i - 4 + 3)%3];
        pprev = buf[(i - 5 + 3)%3];
        //提取出上一行所检测到的角点位置
        cornerpos = cpbuf[(i - 4 + 3)%3];
        //提取出上一行的角点数
        ncorners = cornerpos[-1];
        //在上一行内遍历整个检测到的角点
        for( k = 0; k < ncorners; k++ )
        {
			
            j = cornerpos[k];    //得到角点的位置
            score = prev[j];    //得到该角点的得分函数值
            //在3×3的角点邻域内，计算当前角点是否为最大值，如果是则压入特性值向量中
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
