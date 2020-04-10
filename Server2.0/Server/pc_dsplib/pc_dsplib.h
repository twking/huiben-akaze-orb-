#ifndef _TW_DSPLIB__
#define _TW_DSPLIB__

#include "tw_tpyedef.h"
#include "opencv2/features2d.hpp"
#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/imgproc/types_c.h"



using namespace cv;
#ifdef __cplusplus
extern "C" {
#endif

int initdetectAndDescriptors(INDENTIFOptions *options);
int detectAndDescriptors(cv::Mat img_1, pTEMPLATEDES pkeyAndDes);
void getKeypoint(cv::Mat  image, std::vector<KeyPoint>& tempkeypoints, unsigned int leavekeyNums);
bool SortByM3(const KeyPoint &v1, const KeyPoint &v2);
void viewTest(cv::Mat img_temp,cv::Mat img_src);

#ifdef __cplusplus
}
#endif

#endif