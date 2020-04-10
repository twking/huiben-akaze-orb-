#ifndef _TW_TRAN__
#define _TW_TRAN__

#include "../pc_dsplib/tw_tpyedef.h"
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

int inittranBookDescrib(INDENTIFOptions *options);
int tranBookDescrib(string str, string fileType);
int tranContentDescrib(string str, string fileType);

/*
path: 指定目录
files: 保存结果
fileType: 指定的文件格式，如 .jpg
*/
void getAllFiles(string path, vector<string>& files, string fileType);


#ifdef __cplusplus
}
#endif

#endif