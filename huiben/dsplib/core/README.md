
0、修改交叉编译器
　　打开文件“t20-linux.cmake”，更该“GXX_DIR”
3. 生成Make工程

```
cd build
cmake -DCMAKE_TOOLCHAIN_FILE=../t20-linux.cmake ..
make
```

