#! /bin/bash

mpicc -std=c99 -O3 parallel_sz_2.0.c -o parallel_sz_old -I/home/xin/utils/sz_2.0.0/include /home/xin/utils/sz_2.0.0/lib/libSZ.a /home/xin/utils/sz_2.0.0/lib/libzlib.a -lm

mpicc -std=c99 -O3 parallel_zfp.c -o parallel_zfp -I/home/xin/codes/zfp-0.5.4/include -I/home/xin/utils/sz_master/include /home/xin/utils/sz_master/lib/libSZ.a /home/xin/utils/sz_master/lib/libzlib.a /home/xin/utils/sz_master/lib/libzstd.a /home/xin/codes/zfp-0.5.4/lib/libzfp.a -lm

g++ -O3 -c sz_zfp_selector.cpp -o selector.o -I/home/xin/codes/zfp-0.3.1/inc -I/home/xin/utils/sz_master/include

mpicxx -O3 parallel_selector.c selector.o -o parallel_selector -I/home/xin/codes/zfp-0.3.1/inc -I/home/xin/utils/sz_master/include /home/xin/utils/sz_master/lib/libSZ.a /home/xin/utils/sz_master/lib/libzlib.a /home/xin/utils/sz_master/lib/libzstd.a /home/xin/codes/zfp-0.3.1/lib/libzfp.a

 
