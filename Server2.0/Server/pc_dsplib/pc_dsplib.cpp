#include "pc_dsplib.h"

#define SCR_W 1280
#define SCR_H 720
#define KKK 2
#define DEST_W (565*KKK)
#define DEST_H (346*KKK)
static INDENTIFOptions akaze_options;
int initdetectAndDescriptors(INDENTIFOptions *options)
{
	if (options->dthreshold > 0.0005f&&options->dthreshold < 0.1f)
	{
		akaze_options.dthreshold = options->dthreshold;
	}
	else
	{
		return -1;
	}
	if (options->omax >= 1&& options->omax <= 4)
	{
		akaze_options.omax = options->omax;
	}
	else
	{
		return -1;
	}
	if (options->nsublevels >= 1 && options->nsublevels <= 4)
	{
		akaze_options.nsublevels = options->nsublevels;
	}
	else
	{
		return -1;
	}
	if (options->teamplatekeyPointNums >= 50&& options->teamplatekeyPointNums < 500)
	{
		akaze_options.teamplatekeyPointNums = options->teamplatekeyPointNums;
	}
	else
	{
		return -1;
	}
	if (options->keymode== INDENTIFOptions::DATABASE_KEY_MODE_MAX|| options->keymode == INDENTIFOptions::DATABASE_KEY_MODE_MEAN)
	{
		akaze_options.keymode = options->keymode;
	}
	else
	{
		return -1;
	}
	if (options->keysize>=5 && options->keysize <= 15)
	{
		akaze_options.keysize = options->keysize;
	}
	else
	{
		return -1;
	}

	if (options->img_max_width >= 300 && options->img_max_width <= 1500)
	{
		akaze_options.img_max_width = options->img_max_width;
	}
	else
	{
		return -1;
	}
	if (options->img_max_height >= 300 && options->img_max_height <= 1500)
	{
		akaze_options.img_max_height = options->img_max_height;
	}
	else
	{
		return -1;
	}
	return 1;
}
bool SortByM3(const KeyPoint &v1, const KeyPoint &v2)//注意：本函数的参数的类型一定要与vector中元素的类型一致  
{
	return v1.response > v2.response;//降序排列  
}

void getKeypoint(cv::Mat  image, std::vector<KeyPoint>& tempkeypoints, unsigned int leavekeyNums)
{
	unsigned int i, j, k_size, one_k, n;
	CV_Assert(!image.empty());
	typedef  std::vector<std::vector<std::vector<KeyPoint>>>	PointVector3D;
	PointVector3D *pkeyp3D = new PointVector3D;
	unsigned int w, h;
	h = image.rows / akaze_options.keysize + 1;
	w = image.cols / akaze_options.keysize + 1;
	unsigned int t_w, t_h;
	// 动态设置大小.
	pkeyp3D->resize(h);
	for (i = 0; i<h; ++i)
	{
		(*pkeyp3D)[i].resize(w);
	}

	k_size = tempkeypoints.size();
	for (i = 0; i<tempkeypoints.size(); i++)
	{
		t_w = (unsigned int)(tempkeypoints[i].pt.x / akaze_options.keysize);
		t_h = (unsigned int)(tempkeypoints[i].pt.y / akaze_options.keysize);
		(*pkeyp3D)[t_h][t_w].push_back(tempkeypoints[i]);
	}
	for (i = 0; i<h; ++i)
	{
		for (j = 0; j<w; ++j)
		{
			std::sort((*pkeyp3D)[i][j].begin(), (*pkeyp3D)[i][j].end(), SortByM3);
		}
	}
	one_k = 4;
	while (one_k>1)
	{
		for (i = 0; i<h; ++i)
		{
			for (j = 0; j<w; ++j)
			{
				while ((*pkeyp3D)[i][j].size() >= one_k)
				{
					if (k_size <= leavekeyNums)
					{
						goto RET;
					}
					k_size--;
					(*pkeyp3D)[i][j].pop_back();

				}
			}
		}
		one_k--;
	}

	tempkeypoints.erase(tempkeypoints.begin(), tempkeypoints.end());
	for (i = 0; i<h; ++i)
	{
		for (j = 0; j<w; ++j)
		{
			for (n = 0; n<(*pkeyp3D)[i][j].size(); n++)
			{
				tempkeypoints.push_back((*pkeyp3D)[i][j][n]);
			}
			(*pkeyp3D)[i][j].erase((*pkeyp3D)[i][j].begin(), (*pkeyp3D)[i][j].end());
		}
	}
	k_size = tempkeypoints.size();
	if (k_size>leavekeyNums)
	{
		k_size = tempkeypoints.size();
		for (i = 0; i<tempkeypoints.size(); i++)
		{
			t_w = (unsigned int)((tempkeypoints[i].pt.x + akaze_options.keysize / 2) / akaze_options.keysize);
			t_h = (unsigned int)((tempkeypoints[i].pt.y + akaze_options.keysize / 2) / akaze_options.keysize);
			(*pkeyp3D)[t_h][t_w].push_back(tempkeypoints[i]);
		}
		for (i = 0; i<h; ++i)
		{
			for (j = 0; j<w; ++j)
			{
				std::sort((*pkeyp3D)[i][j].begin(), (*pkeyp3D)[i][j].end(), SortByM3);
			}
		}
		one_k = 4;
		while (one_k>1)
		{
			for (i = 0; i<h; ++i)
			{
				for (j = 0; j<w; ++j)
				{
					while ((*pkeyp3D)[i][j].size() >= one_k)
					{
						if (k_size <= leavekeyNums)
						{
							goto RET;
						}
						k_size--;
						(*pkeyp3D)[i][j].pop_back();

					}
				}
			}
			one_k--;
		}
	}
RET:
	tempkeypoints.erase(tempkeypoints.begin(), tempkeypoints.end());
	for (i = 0; i<h; ++i)
	{
		for (j = 0; j<w; ++j)
		{
			for (n = 0; n<(*pkeyp3D)[i][j].size(); n++)
			{
				tempkeypoints.push_back((*pkeyp3D)[i][j][n]);
			}
		}
	}
	if (k_size >= leavekeyNums)
	{
		std::sort(tempkeypoints.begin(), tempkeypoints.end(), SortByM3);
		tempkeypoints.erase(tempkeypoints.begin() + leavekeyNums, tempkeypoints.end());
	}

	delete pkeyp3D;
}
int detectAndDescriptors(cv::Mat img_2, pTEMPLATEDES pkeyAndDes)
{
	int keynums;
	cv::Mat img_1;
	if (img_2.rows <= 0 || img_2.cols <= 0 )
		return -1;
	if (img_2.cols > akaze_options.img_max_width || img_2.rows > akaze_options.img_max_height)
	{
		Size a;
		float k1, k2;
		k1 = 1.0f*img_2.cols / akaze_options.img_max_width;
		k2 = 1.0f*img_2.rows / akaze_options.img_max_height;
		if (k2 < k1)
			k2 = k1;
		a.height = (int)(img_2.rows / k2);
		a.width = (int)(img_2.cols / k2);
		resize(img_2, img_1, a);
	}
	else
	{
		img_1 = img_2;
	}
	std::vector<KeyPoint> keypoints_dest;
	cv::Mat dstDesc;
#if USING_MODES==T2X_SYSTEM
	Ptr<AKAZE> modes = AKAZE::create(AKAZE::DESCRIPTOR_MLDB, 256, 3, akaze_options.dthreshold, akaze_options.omax, akaze_options.nsublevels, KAZE::DIFF_PM_G2);
#elif USING_MODES==GK_SYSTEM
	Ptr<ORB> modes = ORB::create(3*akaze_options.teamplatekeyPointNums, 1.2f,4, 31, 0, 2, ORB::HARRIS_SCORE, 31,33);
#else 
	#error "USING_MODES must is AKAZE_DEFINE or ORB_DEFINE"
#endif

	modes->detect(img_1, keypoints_dest);
	keynums = (int)keypoints_dest.size();
	if(keynums>= akaze_options.teamplatekeyPointNums)
	{ 
		if (akaze_options.keymode == INDENTIFOptions::DATABASE_KEY_MODE_MAX)
		{
			std::sort(keypoints_dest.begin(), keypoints_dest.end(), SortByM3);
			keypoints_dest.erase(keypoints_dest.begin()+ akaze_options.teamplatekeyPointNums, keypoints_dest.end());
		}
		else
		{
			getKeypoint(img_1, keypoints_dest, akaze_options.teamplatekeyPointNums);
		}
		
		/*
		cv::Mat img_k;
		drawKeypoints(img_1, keypoints_dest, img_k);
		imshow("keypomt", img_k);
		waitKey(0);
		*/
		modes->compute(img_1, keypoints_dest, dstDesc);
		POINYXY *pKey=(POINYXY *)pkeyAndDes->ponit;
		DESCRIP *pDes = (DESCRIP *)pkeyAndDes->descrip;
		for (int j = 0; j<akaze_options.teamplatekeyPointNums; j++)
		{
			pKey[j].x = keypoints_dest[j].pt.x;
			pKey[j].y = keypoints_dest[j].pt.y;
		}
		memcpy(pDes, dstDesc.data, akaze_options.teamplatekeyPointNums * 32);
		return akaze_options.teamplatekeyPointNums;
	}
	else
	{
		modes->compute(img_1, keypoints_dest, dstDesc);
		POINYXY *pKey=(POINYXY *)pkeyAndDes->ponit;
		DESCRIP *pDes = (DESCRIP *)pkeyAndDes->descrip;
		for (int j = 0; j<keynums; j++)
		{
			pKey[j].x = keypoints_dest[j].pt.x;
			pKey[j].y = keypoints_dest[j].pt.y;
		}
		memcpy(pDes, dstDesc.data, keynums * 32);

		printf("warning:ak_detectAndDescriptors keynums = %d<%d\n", keynums,akaze_options.teamplatekeyPointNums);
		return keynums;
	}
	
}


