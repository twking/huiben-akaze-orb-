#ifndef _TW_MKDIR_IO_H__
#define _TW_MKDIR_IO_H__

#ifdef WIN32
#define SD_PIC_DIR "F:/sdcard/pic/"  
#include <io.h>
#include <direct.h>          
#endif  
#ifdef linux
#include <unistd.h>
#include <sys/stat.h> 
#define SD_PIC_DIR "/mnt/sdcard/pic/"        
#endif

#ifdef __cplusplus
extern "C" {
#endif

int CreateDir(const char*sPathName);  


#ifdef __cplusplus
}
#endif

#endif