
#include "matchDestination/pc_huiben.h"
#include "pc_dsplib/pc_dsplib.h"

#pragma warning(disable:4996)

int domeonpagesdesrictofile(string img_str, string filesname)
{

	Mat img_1;
	Mat gimg_1;
	TEMPLATEDES keyAndDes;
	int retcn,i;
	unsigned char *preturnData;
	unsigned int addjiaoyan=0;

	keyAndDes.descrip = new DESCRIP[50];
	keyAndDes.ponit = new POINYXY[50];
	
	FILE *fid;
	fid = fopen(filesname.c_str(), "wb");
	if (fid == NULL)
	{
		printf("write file failed!\n");
		return -1;
	}
	INDENTIFOptions options;
	options.teamplatekeyPointNums = 50;
	options.keymode = INDENTIFOptions::DATABASE_KEY_MODE_MAX;
	options.dthreshold = 0.0007f;
	options.nsublevels = 1;
	options.keysize = 15;
	initdetectAndDescriptors(&options);

	img_1 = cv::imread(img_str);
	if(img_1.data==NULL)
		return -1;
	cvtColor(img_1, gimg_1, cv::COLOR_RGB2GRAY);
	int keynums = detectAndDescriptors(gimg_1, &keyAndDes);
	if(keynums>10)
	{
		retcn = sizeof(int) + 32 * sizeof(char) + keynums * sizeof(DESCRIP) + keynums * sizeof(POINYXY)+sizeof(addjiaoyan);
		preturnData= (unsigned char *)malloc(retcn);
		if (preturnData != NULL)
		{
			pTEMPLATEDES temp = (pTEMPLATEDES)preturnData;			
			memset(preturnData, 0,retcn);			
			temp->pages_key_nums = keynums;
			memcpy(preturnData + sizeof(int) + 32 * sizeof(char), keyAndDes.descrip, keynums * sizeof(DESCRIP));
			memcpy(preturnData + sizeof(int) + 32 * sizeof(char)+keynums * sizeof(DESCRIP), keyAndDes.ponit, keynums * sizeof(POINYXY));
		}
		for(i=0;i<retcn-4;i++)
		{
			addjiaoyan +=preturnData[i];
		}
		preturnData[i++]=(addjiaoyan)&0xff;
		preturnData[i++]=(addjiaoyan>>8)&0xff;
		preturnData[i++]=(addjiaoyan>>16)&0xff;
		preturnData[i++]=(addjiaoyan>>24)&0xff;
		fwrite(preturnData, retcn, 1, fid);
	}
	fclose(fid);
	if (keyAndDes.descrip != NULL) delete[] keyAndDes.descrip;
	if (keyAndDes.ponit != NULL) delete[] keyAndDes.ponit;
	printf("Write file finshed\n");
	return keynums;
}

void match_testmain(void)
{

	string str = "D:\\booklist\\955713de8a6119e1_t";
	BOOKDESCRIP contensDes;
	long tempsize;
	tempsize = readDesripFromfile(str, &contensDes);
	if(tempsize==0) return;
	char *retStr=NULL;
	string img_str= "D:/pic/1.bmp";
	string filesname= "D:/pic/sample_test";

	///////////以下步骤可以等价于从服务器从设备端获取到的的特征数据
	//domeonpagesdesrictofile(img_str, filesname);
	void * pdecvToPCDesribe = NULL;
	FILE *fid;
	fid = fopen(filesname.c_str(), "rb");
	if (fid == NULL)
	{
		printf("read file failed!\n");
		return ;
	}
	//获取文件大小
	fseek(fid, 0, SEEK_END);
	long lSize = ftell(fid);
	rewind(fid);
	//开辟存储空间
	pdecvToPCDesribe = malloc(lSize);
	if (pdecvToPCDesribe == NULL)
	{
		printf("malloc memory spach failed!\n");
		return ;
	}
	fread(pdecvToPCDesribe, sizeof(char), lSize, fid);
	fclose(fid);
	////////////////////////////////////////////////////////


	if (tempsize > 0)
	{
		retStr = pc_identifycontest(pdecvToPCDesribe, &contensDes);
		if (retStr != NULL)
		{
			printf("find：%s\n", retStr);//服务器进行相关操作，上传数据给设备等等。可以通过retStr索引到目标数据库文件和对应的音频文件
		}
	}

	free(pdecvToPCDesribe);
	free(contensDes.pDescrip);

}

void mainM()
{

	match_testmain();
}