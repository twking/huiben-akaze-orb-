#ifndef _TW_IDENTIFYPICTURES__
#define _TW_IDENTIFYPICTURES__

#include "dsplib/tw_tpyedef.h"
//#include "picData.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
���ܣ��޸�identifybooks ��getPagesDescribe���е��ڲ����в�����
	  �ڵ���identifybooks ��getPagesDescribe����ʱ������initIdentify��ʼ��options������
����options:
 INDENTIFOptions()
        : omax(1)
        , nsublevels(1)
        , img_max_width(IMAG_TEMPLATE_W)
        , img_max_height(DEST_TEMPLATE_H)
        , dthreshold(0.0007f)//[0.005f,0.1f]
		, keymode(DATABASE_KEY_MODE_MEAN)
		, keysize(5)//[5,15],һ������������10*10
		, teamplatekeyPointNums(200)//[200,500]��ģ�����ݿ�������������
		, srckeyPointNums(50)//[50,200]�������ļ������������
		, indentCNNums(8)//[8 15] Ҫ����ƥ�䵽������������
		, rectificationflag(1)//ͼ���Ƿ������־��0��ʾ��������1��ʾ����
        , dmatchratio(0.8f)//[0.5,0.85] ��ֵ
		������identifybooks ��getPagesDescribe���е��ڲ�Ĭ�ϲ�����
����ֵ��-1,��ʶ��ʼʧ�ܣ�˵��ĳ���������ò�����Ҫ��1 �ɹ���
*/
int  initIdentify(INDENTIFOptions *options);
/*
���ܣ����鱾�����ļ����鱾�б������ļ���ƥ��pGrayimageָ���ͼƬ�����ƥ��ɹ�������Ӧ�ַ�����ʧ�ܷ���NULL
����
1���� rectificationflag=1ʱ���������ܴ���
pGrayimage->width=SCR_W ;//(picData.h)
pGrayimage->height=SCR_H;//(picData.h)
pGrayimage->imageData =ָ��Ҷ����ݣ�
����û����Ҫ��
2��pBOOKDescribe,�鱾�������ļ��ṹ����ָ��
3��pCONTENTESDescribe���鱾�б�������ļ��ṹ��������
����ֵ��ƥ��ɹ���������Ӧ���ַ���������ַ�������==16��˵������pCONTENTESDescribe�����ҵ�������˵������pBOOKDescribe�����ҵ�
		ƥ��ʧ�ܣ�NULL
		��һ��ƥ��ɹ�����ǰ�κ���һ��ͼƬ��ͬ��SAMEPIC_FIND //"samepic_finded"
		��һ��ƥ��ʧ�ܣ���ǰ�κ���һ��ͼƬ��ͬ���������ƥ��ʧ�ܣ�SAMEPIC_NOTFIND //"samepic_notfind"
*/
const char *identifybooks(Imagetype* pGrayimage, pBOOKDESCRIP pBOOKDescribe, pBOOKDESCRIP pCONTENTESDescribe);
/*
���ܣ���ȡ�����ļ����ɹ����ظ������ļ��ĳ��ȣ�ʧ�ܷ���0��
	 ������سɹ���ע�⵱ʹ����pBooksDesc����ʱ����������²���
	 if (pBooksDesc->pDescrip != NULL)
	 {
		free(pBooksDesc->pDescrip);//�ǵ��ͷ�pBooksDesc->pDescrip�ڴ�ռ�
		pBooksDesc->pDescrip = NULL;
	 }
����
1��str��Ҫ��ȡ�����ļ���Ŀ¼�ַ���
2��pBOOKDescribe,�鱾�������ļ��ṹ����ָ��
3��pCONTENTESDescribe���鱾�б�������ļ��ṹ��������
����ֵ��ƥ��ɹ���������Ӧ���ַ���������ַ�������==16��˵������pCONTENTESDescribe�����ҵ�������˵������pBOOKDescribe�����ҵ�
		ƥ��ʧ�ܣ�NULL
*/
int readDesripFromfile(string str, pBOOKDESCRIP pBooksDesc);
/*
���ܣ�����pGrayimageͼƬ������ֵ
����
1���� rectificationflag=1ʱ���������ܴ���
pGrayimage->width=SCR_W ;//(picData.h)
pGrayimage->height=SCR_H;//(picData.h)
pGrayimage->imageData =ָ��Ҷ����ݣ�
����û����Ҫ��
����ֵ���ɹ���������Ӧ������ָ�룬���ݳ���==*cont
		ʧ�ܣ�NULL
*/
void *getPagesDescribe(Imagetype* pGrayimage,int * cont);
/*
���ܣ���fileBookList�����ļ��������fileDownloadDesc�����ķ����������ݣ�����ڸú��������г����豸ͻȻ�쳣��
	  ��fileBookList���ܱ��𺦡�����֮ǰ�����ñ�����ʩ
����
1��fileBookList��Ҫ���µ��鱾�б������ļ�
2��fileDownloadDesc��Ҫ��ǰ�ķ�����鱾�����ļ�
����ֵ���ɹ������º��fileBookList�ļ���С ��λbytes
		ʧ�ܣ�0
*/
int upBookList(string fileBookList, string fileDownloadDesc);
/*
���ܣ���fileBookList�����ļ�����ɾ���鱾���к�=booksrial_number�ķ�����������
����
1��fileBookList��Ҫ���µ��鱾�б������ļ�
2��booksrial_number���鱾���кż�����������ļ�����16 bytes
����ֵ���ɹ������º��fileBookList�ļ���С ��λbytes
		ʧ�ܣ�0
*/
int delBookList(string fileBookList, string booksrial_number);
/*
���ܣ���fileBookList�����ļ���������鱾���к�=booksrial_number�ķ�����������
����
1��fileBookList��Ҫ���µ��鱾�б������ļ�
2��booksrial_number���鱾���кż�����������ļ�����16 bytes
����ֵ���ɹ���>=1����ʾ���ҵ�booksrial_number���ҳ��ִ���Ϊ����ֵ
		ʧ�ܣ�0
*/
int SearchBooksrial_number(string fileBookList, string booksrial_number);


/*
���ܣ���⵱ǰ��ͼƬ����һ��ͼƬ�Ƿ�����ͬ��
����
1��pImage ָ��Ҷ�ͼƬ����
����ֵ��1����ǰ�κ���һ��ͼƬ��ͬ
		0����ǰ�κ���һ��ͼƬ����ͬ
		-1������ִ��ʧ��
*/
int detetcPIC_three(Imagetype * pImage);

#ifdef __cplusplus
}
#endif

#endif
