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
	inittranBookDescrib(&options);//初始化一些训练参数，按照上面设计即可
	tranBookDescrib(bookstr, fileType);//获取书籍的训练数据接口
	tranContentDescrib(contentstr, contentType);//获取目录册的训练数据接口
}
