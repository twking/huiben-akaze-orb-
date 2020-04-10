#include "pc_huiben.h"

static const uchar popCountTable[256] =
{
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};
inline int twnormHamming_i(const unsigned int* a, const unsigned int* b, int n)
{
    int i = 0;
    int result = 0;
	unsigned int hanm;
	unsigned char *p = (unsigned char *)&hanm;
    for( ; i < n; i++ )
	{
		hanm = a[i] ^ b[i];
        result += popCountTable[p[0]];
		result += popCountTable[p[1]];
		result += popCountTable[p[2]];
		result += popCountTable[p[3]];
	}
    return result;
}
unsigned int hanmmingmatch(pDESCRIP templateDesrcrp,unsigned int templateDesrcrpnums,pDESCRIP srcDesrcrp,unsigned int srcnums,pMatchtype Ptr,float k)
{
	unsigned int i,j,tempdistance,returnval=0;
	Matchtype first,secend;
	first.srcIndex =0;
	first.templateIndex =0;
	first.distance=256;
	secend = first;
	for(i=0;i<srcnums;i++)
	{
		first.srcIndex =i;
		first.templateIndex =0;
		first.distance=256;
		secend = first;
		for(j=0;j<templateDesrcrpnums;j++)
		{
			tempdistance = twnormHamming_i((unsigned int*) templateDesrcrp[j].descrip,(unsigned int*)srcDesrcrp[i].descrip,8);
			if(tempdistance<secend.distance)
			{
				if(tempdistance<first.distance)
				{
					secend = first;
					first.templateIndex = j;
					first.distance = tempdistance;
				}
				else
				{
					secend.templateIndex = j;
					secend.distance = tempdistance;
				}
			}
		}
		if(first.distance<secend.distance*k)
		{
			
			Ptr[returnval] = first;
			returnval ++;
		}
	}
	return returnval;
}
