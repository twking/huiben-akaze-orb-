#ifndef C_BMP_H
#define C_BMP_H

typedef struct
{
    //unsigned short    bfType;   
    unsigned int    bfSize;  
    unsigned short    bfReserved1;
    unsigned short    bfReserved2;
    unsigned int    bfOffBits;  
} ClBitMapFileHeader;  

typedef struct  
{  
    unsigned int  biSize;   
    int   biWidth;   
    int   biHeight;   
    unsigned short   biPlanes;   
    unsigned short   biBitCount;  
    unsigned int  biCompression;   
    unsigned int  biSizeImage;   
    int   biXPelsPerMeter;   
    int   biYPelsPerMeter;   
    unsigned int   biClrUsed;   
    unsigned int   biClrImportant;   
} ClBitMapInfoHeader;  

typedef struct   
{  
    unsigned char rgbBlue; //该颜色的蓝色分量   
    unsigned char rgbGreen; //该颜色的绿色分量   
    unsigned char rgbRed; //该颜色的红色分量   
    unsigned char rgbReserved; //保留值   
} ClRgbQuad;  

typedef struct  
{  
    int width;  
    int height;  
    int channels;  
    unsigned char* imageData;  
}ClImage;  

ClImage* clLoadImage(char* path);  
bool clSaveImage(char* path, ClImage* bmpImg);  

#endif