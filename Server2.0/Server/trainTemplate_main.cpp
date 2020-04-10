#include "trainTemplate/tran.h"

#define BOOKLIST "F:\\VS2012_project\\database\\booklist\\"
#define BOOKS "F:\\VS2012_project\\database\\books\\"
#define TESTPICYURE "F:\\VS2012_project\\database\\testPicture\\"

void main()
{
	INDENTIFOptions options;
	string bookstr = BOOKS;
	bookstr.append("da_wei_re_ma_fang\\Pictures");//wo_ba_ba hai_xiu_de_ai_mi_li da_wei_re_ma_fang
	string fileType = ".jpg";
	string contentType = ".bmp";
	string contentstr = BOOKS;
	contentstr.append("booklist\\Pictures");
	options.dthreshold = 0.001f;
	options.nsublevels = 2;
	inittranBookDescrib(&options);//��ʼ��һЩѵ������������������Ƽ���
	tranBookDescrib(bookstr, fileType);//��ȡ�鼮��ѵ�����ݽӿ�
	tranContentDescrib(contentstr, contentType);//��ȡĿ¼���ѵ�����ݽӿ�
}
