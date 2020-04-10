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
path: ָ��Ŀ¼
files: ������
fileType: ָ�����ļ���ʽ���� .jpg
*/
void getAllFiles(string path, vector<string>& files, string fileType);


#ifdef __cplusplus
}
#endif

#endif