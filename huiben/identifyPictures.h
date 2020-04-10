#ifndef _TW_IDENTIFYPICTURES__
#define _TW_IDENTIFYPICTURES__

#include "dsplib/tw_tpyedef.h"
//#include "picData.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
功能：修改identifybooks ，getPagesDescribe运行的内部运行参数。
	  在调用identifybooks ，getPagesDescribe运行时，先用initIdentify初始化options变量。
参数options:
 INDENTIFOptions()
        : omax(1)
        , nsublevels(1)
        , img_max_width(IMAG_TEMPLATE_W)
        , img_max_height(DEST_TEMPLATE_H)
        , dthreshold(0.0007f)//[0.005f,0.1f]
		, keymode(DATABASE_KEY_MODE_MEAN)
		, keysize(5)//[5,15],一个特征点的面积10*10
		, teamplatekeyPointNums(200)//[200,500]，模板数据库的特征点的数量
		, srckeyPointNums(50)//[50,200]，查找文件特征点的数量
		, indentCNNums(8)//[8 15] 要求能匹配到的特征点数量
		, rectificationflag(1)//图像是否矫正标志，0表示不矫正，1表示矫正
        , dmatchratio(0.8f)//[0.5,0.85] 比值
		以上是identifybooks ，getPagesDescribe运行的内部默认参数。
返回值：-1,标识初始失败，说明某个变量设置不符合要求，1 成功。
*/
int  initIdentify(INDENTIFOptions *options);
/*
功能：从书本特征文件和书本列表特征文件中匹配pGrayimage指向的图片，如果匹配成功返回响应字符串，失败返回NULL
参数
1：当 rectificationflag=1时，矫正功能打开则：
pGrayimage->width=SCR_W ;//(picData.h)
pGrayimage->height=SCR_H;//(picData.h)
pGrayimage->imageData =指向灰度数据，
否则没有这要求
2：pBOOKDescribe,书本的特征文件结构数据指针
3：pCONTENTESDescribe，书本列表的特征文件结构数据制作
返回值：匹配成功：返回响应的字符串，如果字符串长度==16，说明是在pCONTENTESDescribe里面找到，否则说明是在pBOOKDescribe里面找到
		匹配失败：NULL
		上一次匹配成功，当前次和上一次图片相同：SAMEPIC_FIND //"samepic_finded"
		上一次匹配失败，当前次和上一次图片相同，并且这次匹配失败：SAMEPIC_NOTFIND //"samepic_notfind"
*/
const char *identifybooks(Imagetype* pGrayimage, pBOOKDESCRIP pBOOKDescribe, pBOOKDESCRIP pCONTENTESDescribe);
/*
功能：读取特征文件，成功返回该特征文件的长度，失败返回0。
	 如果返回成功则，注意当使用完pBooksDesc变量时，请进行如下操作
	 if (pBooksDesc->pDescrip != NULL)
	 {
		free(pBooksDesc->pDescrip);//记得释放pBooksDesc->pDescrip内存空间
		pBooksDesc->pDescrip = NULL;
	 }
参数
1：str：要读取特征文件的目录字符串
2：pBOOKDescribe,书本的特征文件结构数据指针
3：pCONTENTESDescribe，书本列表的特征文件结构数据制作
返回值：匹配成功：返回响应的字符串，如果字符串长度==16，说明是在pCONTENTESDescribe里面找到，否则说明是在pBOOKDescribe里面找到
		匹配失败：NULL
*/
int readDesripFromfile(string str, pBOOKDESCRIP pBooksDesc);
/*
功能：计算pGrayimage图片的特征值
参数
1：当 rectificationflag=1时，矫正功能打开则：
pGrayimage->width=SCR_W ;//(picData.h)
pGrayimage->height=SCR_H;//(picData.h)
pGrayimage->imageData =指向灰度数据，
否则没有这要求
返回值：成功：返回响应的数据指针，数据长度==*cont
		失败：NULL
*/
void *getPagesDescribe(Imagetype* pGrayimage,int * cont);
/*
功能：向fileBookList特征文件里面添加fileDownloadDesc包含的封面特征数据，如果在该函数运行中出现设备突然异常，
	  则fileBookList可能被损害。调用之前请做好保护措施
参数
1：fileBookList：要更新的书本列表特征文件
2：fileDownloadDesc：要提前的封面的书本特征文件
返回值：成功：更新后的fileBookList文件大小 单位bytes
		失败：0
*/
int upBookList(string fileBookList, string fileDownloadDesc);
/*
功能：向fileBookList特征文件里面删除书本序列号=booksrial_number的封面特征数据
参数
1：fileBookList：要更新的书本列表特征文件
2：booksrial_number：书本序列号即该书的特征文件名称16 bytes
返回值：成功：更新后的fileBookList文件大小 单位bytes
		失败：0
*/
int delBookList(string fileBookList, string booksrial_number);
/*
功能：向fileBookList特征文件里面查找书本序列号=booksrial_number的封面特征数据
参数
1：fileBookList：要更新的书本列表特征文件
2：booksrial_number：书本序列号即该书的特征文件名称16 bytes
返回值：成功：>=1，表示查找到booksrial_number，且出现次数为返回值
		失败：0
*/
int SearchBooksrial_number(string fileBookList, string booksrial_number);


/*
功能：检测当前次图片和上一次图片是否是相同的
参数
1：pImage 指向灰度图片数据
返回值：1：当前次和上一次图片相同
		0：当前次和上一次图片不相同
		-1：函数执行失败
*/
int detetcPIC_three(Imagetype * pImage);

#ifdef __cplusplus
}
#endif

#endif
