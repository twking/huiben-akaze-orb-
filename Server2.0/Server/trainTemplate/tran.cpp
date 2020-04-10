
#include "../pc_dsplib/pc_dsplib.h"
#include "tran.h"
#include "xxhash.h"
#include <ctime>
#include <string.h>
#include <fstream>
#include <io.h>
#include <stdio.h>

#pragma warning(disable:4996)

using namespace std;
using namespace cv;


static INDENTIFOptions tran_options;
static unsigned long long calcul_hash(const void *buffer, size_t length)
{
	unsigned long long const seed = 0; /* or any other value */
	unsigned long long const hash = XXH64(buffer, length, seed);
	return hash;
}
void getAllFiles(string path, vector<string>& files, string fileType)
{
	// 文件句柄
	long hFile = 0;
	// 文件信息
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*"+ fileType).c_str(), &fileinfo)) != -1) {
		do {
			// 保存文件的全路径
			files.push_back(p.assign(path).append("\\").append(fileinfo.name));

		} while (_findnext(hFile, &fileinfo) == 0); //寻找下一个，成功返回0，否则-1

		_findclose(hFile);
	}
}
int inittranBookDescrib(INDENTIFOptions *options)
{
	if (options->dthreshold > 0.0005f&&options->dthreshold < 0.1f)
	{
		tran_options.dthreshold = options->dthreshold;
	}
	else
	{
		return -1;
	}
	if (options->omax >= 1 && options->omax <= 4)
	{
		tran_options.omax = options->omax;
	}
	else
	{
		return -1;
	}
	if (options->nsublevels >= 1 && options->nsublevels <= 4)
	{
		tran_options.nsublevels = options->nsublevels;
	}
	else
	{
		return -1;
	}
	if (options->teamplatekeyPointNums >= 200 && options->teamplatekeyPointNums <= 500)
	{
		tran_options.teamplatekeyPointNums = options->teamplatekeyPointNums;
	}
	else
	{
		return -1;
	}
	if (options->keymode == INDENTIFOptions::DATABASE_KEY_MODE_MAX || options->keymode == INDENTIFOptions::DATABASE_KEY_MODE_MEAN)
	{
		tran_options.keymode = options->keymode;
	}
	else
	{
		return -1;
	}
	if (options->keysize >= 5 && options->keysize <= 15)
	{
		tran_options.keysize = options->keysize;
	}
	else
	{
		return -1;
	}

	if (options->img_max_width >= 300 && options->img_max_width <= 1500)
	{
		tran_options.img_max_width = options->img_max_width;
	}
	else
	{
		return -1;
	}
	if (options->img_max_height >= 300 && options->img_max_height <= 1500)
	{
		tran_options.img_max_height = options->img_max_height;
	}
	else
	{
		return -1;
	}
	return initdetectAndDescriptors(&tran_options);

}
int tranBookDescrib(string str, string fileType)
{
	vector<string> files;
	vector<string> res;
	getAllFiles(str, files, fileType);
	if (files.size() < 1)
	{
		return -1;
	}
	const char * split = "\\";
	char * p;
	//先将要切割的字符串从string类型转换为char*类型  
	char * strs = new char[files[0].length() + 1]; //不要忘了  
	if (strs == NULL)
	{
		return -1;
	}
	strcpy(strs, files[0].c_str());
	p = strtok(strs, split);
	while (p != NULL)
	{
		string s = p; //分割得到的字符串转换为string类型  
		res.push_back(s); //存入结果数组  
		p = strtok(NULL, split);
	}
	delete[] strs;
	if (res.size() <3)
	{
		return -1;
	}
	strs = new char[res[res.size() - 3].length() + 1];
	if (strs == NULL)
	{
		return -1;
	}
	strcpy(strs, res[res.size() - 3].c_str());
	BOOKDESCRIP booksDesr;
	booksDesr.filesize = 0;
	booksDesr.nums_desc = files.size();

	char has128[32] = {0};
	{
		unsigned long long hash_value = calcul_hash((void *)strs, strlen(strs));
		_snprintf(has128, 32, "%llx", hash_value);
		printf("BOOKNAME = %s,hash_value=%lld,has128=%s\n", strs, hash_value, has128);
	}
	memcpy((void *)booksDesr.srial_number, has128, 32);

	delete[] strs;
	string tempstr,istr;
	Mat img_1;
	Mat gimg_1;
	TEMPLATEDES keyAndDes;
	keyAndDes.descrip = new DESCRIP[tran_options.teamplatekeyPointNums];
	keyAndDes.ponit = new POINYXY[tran_options.teamplatekeyPointNums];
	if (keyAndDes.descrip == NULL|| keyAndDes.ponit ==NULL) return -1;
	booksDesr.pDescrip = &keyAndDes;

	string filestr = str + "\\"+has128;
	FILE *fid;
	fid = fopen(filestr.c_str(), "wb");
	if (fid == NULL)
	{
		printf("write file failed!\n");
		return -1;
	}
	printf("Writing %s file\n", filestr.c_str());
	fwrite((char *)&booksDesr, sizeof(booksDesr.versionType)+sizeof(booksDesr.filesize)+sizeof(booksDesr.nums_desc)+sizeof(booksDesr.srial_number), 1, fid);

	unsigned int lSize = ftell(fid);
	for (unsigned int i = 1; i <= files.size(); i++)
	{
		stringstream num_str;
		num_str << i;
		istr = num_str.str();
		tempstr = str +"\\"+istr + fileType;
		printf("running%s\n",tempstr.c_str());
		img_1 = cv::imread(tempstr);
		cvtColor(img_1, gimg_1, cv::COLOR_RGB2GRAY);
		if (i == 1)
		{
			filestr += ".bmp";
			cv::imwrite(filestr, img_1);
		}
		
		//gimg_1 = img_1;
		int keynums = detectAndDescriptors(gimg_1, &keyAndDes);
		int  index = tempstr.find_last_of('.');
		int start  = tempstr.find_last_of('\\');
		string substring = tempstr.substr(start+1, index - start-1);
		memset(keyAndDes.pages_srial_number,0,32);
		memcpy(keyAndDes.pages_srial_number, substring.c_str(), substring.length());
		keyAndDes.pages_key_nums = keynums;
		fwrite((char *)&keyAndDes, sizeof(keyAndDes.pages_key_nums)+sizeof(keyAndDes.pages_srial_number), 1, fid);
		fwrite(keyAndDes.descrip, sizeof(DESCRIP), keynums, fid);
		fwrite(keyAndDes.ponit, sizeof(POINYXY), keynums, fid);
	}
	unsigned int lastlSize = ftell(fid);
	booksDesr.filesize = lastlSize - lSize;
	rewind(fid);
	fwrite((char *)&booksDesr,sizeof(booksDesr.versionType)+sizeof(booksDesr.filesize), 1, fid);
	fclose(fid);
	if (keyAndDes.descrip != NULL) delete[] keyAndDes.descrip;
	if (keyAndDes.ponit != NULL) delete[] keyAndDes.ponit;
	printf("Write file finshed\n");
	return files.size();


	
}
int tranContentDescrib(string str, string fileType)
{
	vector<string> files;
	vector<string> res;
	getAllFiles(str, files, fileType);
	if (files.size() < 1)
	{
		return -1;
	}
	
	BOOKDESCRIP booksDesr;
	booksDesr.filesize = 0;
	booksDesr.nums_desc = files.size();
	const char *bookuuid = "booklist";
	char has128[32] = {0};
	{
		unsigned long long hash_value = calcul_hash((void *)bookuuid, strlen(bookuuid));
		_snprintf(has128, 32, "%llx", hash_value);
		printf("booklist = %s,hash_value=%lld,has128=%s\n", bookuuid, hash_value, has128);
	}
	memcpy((void *)booksDesr.srial_number, has128, 32);
	Mat img_1;
	Mat gimg_1;
	TEMPLATEDES keyAndDes;
	keyAndDes.descrip = new DESCRIP[tran_options.teamplatekeyPointNums];
	keyAndDes.ponit = new POINYXY[tran_options.teamplatekeyPointNums];
	if (keyAndDes.descrip == NULL || keyAndDes.ponit == NULL) return -1;
	booksDesr.pDescrip = &keyAndDes;

	string filestr = str + "\\" + has128;
	FILE *fid;
	fid = fopen(filestr.c_str(), "wb");
	if (fid == NULL)
	{
		printf("write file failed!\n");
		return -1;
	}
	printf("Writing %s file\n", filestr.c_str());
	fwrite((char *)&booksDesr, sizeof(booksDesr.versionType)+sizeof(booksDesr.filesize)+sizeof(booksDesr.nums_desc)+sizeof(booksDesr.srial_number), 1, fid);
	unsigned int lSize = ftell(fid);

	for (unsigned int i = 0; i < files.size(); i++)
	{
	
		printf("running%s\n",files[i].c_str());
		img_1 = cv::imread(files[i]);
		cvtColor(img_1, gimg_1, cv::COLOR_RGB2GRAY);
		//gimg_1 = img_1;
		int keynums = detectAndDescriptors(gimg_1, &keyAndDes);
		int  index = files[i].find_last_of('.');
		int start = files[i].find_last_of('\\');
		string substring = files[i].substr(start + 1, index - start - 1);
		memset(keyAndDes.pages_srial_number, 0, 32);
		memcpy(keyAndDes.pages_srial_number, substring.c_str(), substring.length());
		keyAndDes.pages_key_nums = keynums;
		fwrite((char *)&keyAndDes, sizeof(keyAndDes.pages_key_nums) + sizeof(keyAndDes.pages_srial_number), 1, fid);
		fwrite(keyAndDes.descrip, sizeof(DESCRIP), keynums, fid);
		fwrite(keyAndDes.ponit, sizeof(POINYXY), keynums, fid);
	}
	unsigned int lastlSize = ftell(fid);
	booksDesr.filesize = lastlSize - lSize;
	rewind(fid);
	fwrite((char *)&booksDesr,sizeof(booksDesr.versionType)+sizeof(booksDesr.filesize), 1, fid);
	fclose(fid);
	if (keyAndDes.descrip != NULL) delete[] keyAndDes.descrip;
	if (keyAndDes.ponit != NULL) delete[] keyAndDes.ponit;
	printf("Write file finshed\n");
	return files.size();
}



