#ifndef _sz_zfp_selector_h
#define _sz_zfp_selector_h

#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstddef>
#include "sz.h"
#include "rw.h"
#include "zfp.h"
#include "bitstream.h"
#include <iostream>

#define MAX(a, b) a>b?a:b
#define MIN(a, b) a<b?a:b
#define BLOCK_SIZE 4

unsigned char * compress_block(float * data, size_t r3, size_t r2, size_t r1, double eb, size_t * out_size, size_t * out_size_before_lossless, int * select);
float * decompress_block(unsigned char * comp_data, size_t comp_data_size, size_t comp_data_size_before_lossless, int select, int r3, int r2, int r1);

#endif