#include "mkdir_io.h"
#include <string.h>
#include <stdio.h>

int CreateDir(const char*sPathName)
{

	int i,len;
	char   DirName[256];  
	len=strlen(sPathName);
	if(len>254)
	{
		printf("CreateDir:strlen(sPathName)>254%d\n",len);
		return -1;
	}
	strcpy(DirName,sPathName);  

	if(DirName[len-1]!='/')  
	strcat(DirName,   "/");  
	len = strlen(DirName);  
	for(i=1;   i<len;   i++)  
	{  
		if(DirName[i]=='/')  
		{  
			DirName[i]   =   0;  
			if(  access(DirName, 0)!=0)  
			{  
			#ifdef WIN32  
				if(mkdir(DirName)==-1)
				{

					perror("mkdir   error");   
					return   -1;   
				} 
			#endif  
			#ifdef linux   
				if(mkdir(DirName, 0777)==-1)
				{

					perror("mkdir   error");   
					return   -1;   
				}
			#endif 
  
			}  
			DirName[i]   =   '/';  
		}  
	}  
   
  return   0;  

}