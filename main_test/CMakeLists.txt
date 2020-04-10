
# 最低版本要求，必须
#cmake_minimum_required(VERSION 2.8)

# 工程名
#project(main_test)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wstrict-prototypes -fwrapv -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lm")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ldl")
# 添加include目录,本例中文件在./下,所以直接取${PROJECT_SOURCE_DIR}/
include_directories(${PROJECT_SOURCE_DIR}/include)
#include_directories(${PROJECT_SOURCE_DIR}/img/include)
#include_directories(${PROJECT_SOURCE_DIR}/bmp/include)
#include_directories(${PROJECT_SOURCE_DIR}/dsp/include)
#include_directories(${PROJECT_SOURCE_DIR}/bmp/include)
#include_directories(${PROJECT_SOURCE_DIR}/core/include)
#include_directories(${PROJECT_SOURCE_DIR}/img/include)
# 设置生成的库文件目录,这里设置成和库目录相同
#set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
#设置可执行文件目录,同样设成和库目录相同
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
# 设置C编译器编译选项
set(CMAKE_C_FLAGS "-O2 -DNDEBUG -Wall -pipe -fmessage-length=0 -ffunction-sections -fdata-sections")
#set(CMAKE_C_FLAGS_DEBUG "-g -ggdb -DDEBUG")
#set(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")

#设置C编译器 连接选项
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-gc-sections")

#add_subdirectory(Akaze)
#add_subdirectory(core)
#add_subdirectory(img)
#add_subdirectory(dsp)

link_directories(${PROJECT_BINARY_DIR}/lib)

AUX_SOURCE_DIRECTORY(. DIR_SRCS)
ADD_EXECUTABLE(main ${DIR_SRCS})
#下面的 libxxx 是 第三方库 的名字
target_link_libraries(main libdsp.a)
#target_link_libraries(main_test libbmp.a)
#target_link_libraries(main_test libimg.a)
#target_link_libraries(main_test libcore.a)

# 添加子目录（每个子目录下都必须有自已的CMakeLists.txt文件)


