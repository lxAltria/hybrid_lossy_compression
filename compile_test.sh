#! /bin/bash

g++ -O3 -c sz_zfp_selector.cpp -o selector.o -I/home/xin/codes/zfp-0.3.1/inc -I/home/xin/utils/sz_master/include

g++ -O3 sz_zfp_select.cpp selector.o -o select -I/home/xin/codes/zfp-0.3.1/inc -I/home/xin/utils/sz_master/include /home/xin/utils/sz_master/lib/libSZ.a /home/xin/utils/sz_master/lib/libzlib.a /home/xin/utils/sz_master/lib/libzstd.a /home/xin/codes/zfp-0.3.1/lib/libzfp.a
