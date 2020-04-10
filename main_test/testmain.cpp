
#include <ctime>
#include "testmain.h"
#include <string.h>
#pragma warning(disable:4996)
using namespace std;
clock_t start, finish; 
double   duration; 
void appidentifybooks_Books(void)
{
	printf("\n*********************************\n");
	printf("\appidentifybooks_Books test:\n");
	string str = "D:\\booklist\\121f110cbf98ac5b";
	char picfiles[128] = { 0 };
	const char* pages;
	unsigned char *pGrayimage;
	BOOKDESCRIP booksDes;
	booksDes.pDescrip = NULL;
	BOOKDESCRIP contensDes;
	contensDes.pDescrip = NULL;
	long tempsize;
	INDENTIFOptions options;
	initIdentify(&options);
	tempsize = readDesripFromfile(str,&contensDes);
	
	if (tempsize != NULL)
	{
		printf("contensDes->filesize=%d\n", contensDes.filesize);
		printf("contensDes->nums_desc=%d\n", contensDes.nums_desc);
		char tstr[33]={0};
		memset(tstr,0,33);
		printf("strlen(contensDes.srial_number)=%d\n",strlen(contensDes.srial_number));
		memcpy(tstr,contensDes.srial_number, 32);
		printf("contensDes->srial_number=%s\n", tstr);
	}
	else
	{
		printf("获取书籍目录的特征描述失败!\n");
		return;
	}
	int i = 1;
	while (1)
	{
		sprintf(picfiles, "D:/pic/%d.bmp", i);//"D:/pic/y_%d.bmp"
		printf("searching p1_bmp:%s\n",picfiles);
		i++;
		ClImage* img = clLoadImage(picfiles);
		if(i>500)
		{
			break;
		}
		if(img==NULL)
		{
			printf("pic:%s is nor find\n",picfiles);
			continue;
		}

		Imagetype srcImag;
		srcImag.height=img->height;
		srcImag.width=img->width;
		srcImag.imageData=img->imageData;
		if (srcImag.imageData == NULL)
		{
			printf("要识别的图片：%s 获取失败！\n", picfiles);
			if(booksDes.pDescrip!=NULL)
			{ 
				free(booksDes.pDescrip);//记得释放pbooksDes内存空间
				booksDes.pDescrip = NULL;
			}
			if (contensDes.pDescrip != NULL)
			{
				free(contensDes.pDescrip);//记得释放pbooksDes内存空间
				contensDes.pDescrip = NULL;
			}
			return;
		}
		//以上代码模拟摄像头获取到的灰度数据pGrayimage
		pages = identifybooks(&srcImag, &booksDes, &contensDes);//identifybooks(&srcImag,&booksDes, &contensDes);
		//pages =NULL;
		if (pages != NULL)//pages>=0时说明找到匹配图像
		{

			if (strlen(pages) == 16)//表示是在目录数据库找到的
			{
				if (booksDes.pDescrip != NULL)
				{
					free(booksDes.pDescrip);//记得释放pbooksDes内存空间
					booksDes.pDescrip = NULL;
				}
				tempsize = readDesripFromfile(pages,&booksDes);
				if (tempsize != 0)
				{
					printf("get books description data\n");
					printf("contensDes->filesize=%d\n", booksDes.filesize);
					char tstr[33];
					memset(tstr,0,33);
					memcpy(tstr,booksDes.srial_number, 32);
					printf("contensDes->srial_number=%s\n", tstr);

				}
			}
			else//表示是在书籍数据库找到的
			{
				//可以进行播放相关音乐
				char tstr[33];
				memset(tstr,0,33);
				memcpy(tstr, booksDes.srial_number, 32);
				printf("find pic in %s pages=%s\n", tstr, pages);
				//printf("find pic pages=%s\n", pages);
			}
		}
		else
		{
			printf("not find pic:%s\n",picfiles);
		}
		if(img->imageData!=NULL)
		{
			//free(img->imageData);
			img->imageData = NULL;
		}
		if(img!=NULL)
		{
			free(img);
			img = NULL;
		}
	}
	if (booksDes.pDescrip != NULL)
	{
		free(booksDes.pDescrip);//记得释放pbooksDes内存空间
		booksDes.pDescrip = NULL;
	}
	if (contensDes.pDescrip != NULL)
	{
		free(contensDes.pDescrip);//记得释放pbooksDes内存空间
		contensDes.pDescrip = NULL;
	}

}
void appgetPagesDescribe(void)
{
	char *pretval;
	int cn=0;
	char *img_str = "p1_bmp\\y_43.bmp";
	string filesname = "sample_test";
	ClImage* img = clLoadImage(img_str);
	if (img == NULL)
	{
		printf("pic:%s is nor find\n", img_str);
	}

	Imagetype srcImag;
	srcImag.height = img->height;
	srcImag.width = img->width;
	srcImag.imageData = img->imageData;
	pretval = (char *)getPagesDescribe(&srcImag, &cn);

	//
	FILE *fid;
	fid = fopen(filesname.c_str(), "wb");
	if (fid == NULL)
	{
		printf("write file failed!\n");
		return ;
	}
	fwrite(pretval, 1, cn, fid);
	fclose(fid);

	if(pretval!=NULL)//用完记得释放pretval
	{
		free(pretval);
		pretval = NULL;
	}
}
void app_upBookList(void)
{
	string fileBookList="F:\\libdsp_mips_1 .2\\vs2012\\main_test\\booklist\\121f110cbf98ac5b";
	string fileDownloadDesc="F:\\libdsp_mips_1 .2\\vs2012\\main_test\\booklist\\a0b09ca5583c555b";
	long stb = upBookList(fileBookList,fileDownloadDesc);
	printf("app_upBookList :%ld\n",stb);
}
void app_delBookList(void)
{
	string fileBookList="F:\\libdsp_mips_1 .2\\vs2012\\main_test\\booklist\\121f110cbf98ac5b";
	string booksrial_number="a0b09ca5583c555b";
	
	long stb = delBookList(fileBookList,booksrial_number);
	printf("app_delBookList :%ld\n",stb);

}
void appdetetcPIC_three(void)
{
	printf("\n*********************************\n");
	printf("\detetcPIC_three test:\n");

	char picfiles[128] = { 0 };
	int cn;
	unsigned char *pGrayimage;
	int i = 1;
	while (1)
	{
		sprintf(picfiles, "pic/y_%d.bmp", i);
		printf("searching pic:%s\n",picfiles);
		i++;
		ClImage* img = clLoadImage(picfiles);
		if(i>500)
		{
			break;
		}
		if(img==NULL)
		{
			printf("pic:%s is nor find\n",picfiles);
			continue;
		}

		Imagetype srcImag;
		srcImag.height=img->height;
		srcImag.width=img->width;
		srcImag.imageData=img->imageData;
		if (srcImag.imageData == NULL)
		{
			printf("要识别的图片：%s 获取失败！\n", picfiles);
			return;
		}
		//以上代码模拟摄像头获取到的灰度数据pGrayimage
		start = clock(); 
		cn = detetcPIC_three(&srcImag);//identifybooks(&srcImag,&booksDes, &contensDes);
		finish = clock(); 
		duration = (double)(finish - start)/CLOCKS_PER_SEC*1000; 
		
		printf( "detetcPIC_three %f ms\n", duration ); 

		if(img->imageData!=NULL)
		{
			free(img->imageData);
			img->imageData = NULL;
		}
		if(img!=NULL)
		{
			free(img);
			img = NULL;
		}
	}

}
int main()
{
	appidentifybooks_Books();
	//appdetetcPIC_three();
	//CreateDir(SD_PIC_DIR);
	//app_upBookList();
	//app_delBookList();
	//appgetPagesDescribe();
	return 1;
}

