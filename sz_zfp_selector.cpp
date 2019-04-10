#include "sz_zfp_selector.h"

int block_id[4][4][4];

void init_blockid(){
	block_id[0][0][0] = 0;
	block_id[1][0][0] = 1;
	block_id[0][1][0] = 2;
	block_id[0][0][1] = 3;
	block_id[0][1][1] = 4;
	block_id[1][0][1] = 5;
	block_id[1][1][0] = 6;
	block_id[2][0][0] = 7;
	block_id[0][2][0] = 8;
	block_id[0][0][2] = 9;
	block_id[1][1][1] = 10;
	block_id[2][1][0] = 11;
	block_id[2][0][1] = 12;
	block_id[0][2][1] = 13;
	block_id[1][2][0] = 14;
	block_id[1][0][2] = 15;
	block_id[0][1][2] = 16;
	block_id[3][0][0] = 17;
	block_id[0][3][0] = 18;
	block_id[0][0][3] = 19;
	block_id[2][1][1] = 20;
	block_id[1][2][1] = 21;
	block_id[1][1][2] = 22;
	block_id[0][2][2] = 23;
	block_id[2][0][2] = 24;
	block_id[2][2][0] = 25;
	block_id[3][1][0] = 26;
	block_id[3][0][1] = 27;
	block_id[0][3][1] = 28;
	block_id[1][3][0] = 29;
	block_id[1][0][3] = 30;
	block_id[0][1][3] = 31;
	block_id[1][2][2] = 32;
	block_id[2][1][2] = 33;
	block_id[2][2][1] = 34;
	block_id[3][1][1] = 35;
	block_id[1][3][1] = 36;
	block_id[1][1][3] = 37;
	block_id[3][2][0] = 38;
	block_id[3][0][2] = 39;
	block_id[0][3][2] = 40;
	block_id[2][3][0] = 41;
	block_id[2][0][3] = 42;
	block_id[0][2][3] = 43;
	block_id[2][2][2] = 44;
	block_id[3][2][1] = 45;
	block_id[3][1][2] = 46;
	block_id[1][3][2] = 47;
	block_id[2][3][1] = 48;
	block_id[2][1][3] = 49;
	block_id[1][2][3] = 50;
	block_id[0][3][3] = 51;
	block_id[3][0][3] = 52;
	block_id[3][3][0] = 53;
	block_id[3][2][2] = 54;
	block_id[2][3][2] = 55;
	block_id[2][2][3] = 56;
	block_id[1][3][3] = 57;
	block_id[3][1][3] = 58;
	block_id[3][3][1] = 59;
	block_id[2][3][3] = 60;
	block_id[3][2][3] = 61;
	block_id[3][3][2] = 62;
	block_id[3][3][3] = 63;
}

void fp_fwd_lift(signed int* p, unsigned int s)
{
  signed int x = *p; p += s;
  signed int y = *p; p += s;
  signed int z = *p; p += s;
  signed int w = *p; p += s;

  // default, non-orthogonal transform (preferred due to speed and quality)
  //        ( 4  4  4  4) (x)
  // 1/16 * ( 5  1 -1 -5) (y)
  //        (-4  4  4 -4) (z)
  //        (-2  6 -6  2) (w)
  x += w; x >>= 1; w -= x; 
  z += y; z >>= 1; y -= z; 
  x += z; x >>= 1; z -= x; 
  w += y; w >>= 1; y -= w;
  w += y >> 1; y -= w >> 1;

  p -= s; *p = w;
  p -= s; *p = z;
  p -= s; *p = y;
  p -= s; *p = x;
}

void fp_inv_lift(signed int* p, unsigned int s)
{
  signed int x = *p; p += s;
  signed int y = *p; p += s;
  signed int z = *p; p += s;
  signed int w = *p; p += s;

  // default, non-orthogonal transform (preferred due to speed and quality)
  //       ( 4  6 -4 -1) (x)
  // 1/4 * ( 4  2  4  5) (y)
  //       ( 4 -2  4 -5) (z)
  //       ( 4 -6 -4  1) (w)
  y += w >> 1; w -= y >> 1;
  y += w; w <<= 1; w -= y;
  z += x; x <<= 1; x -= z;
  y += z; z <<= 1; z -= y;
  w += x; x <<= 1; x -= w;

  p -= s; *p = w;
  p -= s; *p = z;
  p -= s; *p = y;
  p -= s; *p = x;
}

#if defined(__GNUC__)
#elif defined(__IBMCPP__)
  #include <builtins.h>
#elif defined(_WIN64)
  #include <intrin.h>
  #ifndef HAVE_C99_MATH
    // for old versions of MSVC that do not have C99 math support
    inline long int lrint(double x) { return  (long int)x; }
    inline long long int llrint(double x) { return (long long int)x; }
  #endif
#else
  #error "compiler not supported"
#endif

inline uint fp_uclz(uint32 x)
{
#if defined(__GNUC__)
  return __builtin_clz(x);
#elif defined(__IBMCPP__)
  return __cntlz4(x);
#elif defined(_WIN64)
  unsigned long n;
  _BitScanReverse(&n, x);
  return 31 - n;
#endif
}

inline uint fp_ufls(uint32 x){
#if defined(__GNUC__) || defined(_WIN64)
  return x ? CHAR_BIT * sizeof(x) - fp_uclz(x) : 0;
#elif defined(__IBMCPP__)
  return CHAR_BIT * sizeof(x) - fp_uclz(x);
#endif
}

static uint fp_width(const signed int* data, uint n, uint& m)
{
	while (n--) {
	  signed int x = *data++;
	  m |= x < 0 ? -uint(x) : +uint(x);
	}
	return fp_ufls(m);
}

inline void embedded_encoding(MemoryBitStream& bitstream, const signed int * data, uint maxprec, uint64 count, uint64 width)
{
  MemoryBitStream stream = bitstream;
  uint intprec = 32;
  uint kmin = intprec > maxprec ? intprec - maxprec : 0;
  // output one bit plane at a time from MSB to LSB
  for (uint k = intprec, n = 0; k-- > kmin;) {
    // encode bit plane k
    for (uint i = 0;;) {
      // encode group of significant values
      for (; i < n; i++) {
        // encode bit k of data[i]
        bool sign = data[i] < 0;
        uint x = uint(sign ? -data[i] : +data[i]) >> k;
        // write bit k of |data[i]|
        stream.write(x & uint(1));
        if (x == uint(1)) {
          // write sign bit also
          stream.write(sign);
          // num_sign_bits ++;
        }
      }
      // have all groups been encoded?
      if (!count)
        break;
      // test next group
      if ((width & 0x3fu) > k) {
        // group is significant; peel off and encode first subgroup
        stream.write(true);
        n += count & 0xfu;
        count >>= 4;
        width >>= 6;
      }
      else {
        // group is insignificant; continue with next bit plane
        stream.write(false);
        break;
      }
    }
  }

exit:
  bitstream = stream;
}

inline void embedded_decoding(MemoryBitStream& bitstream, signed int* data, uint maxprec, uint64 count, uint size)
{
  MemoryBitStream stream = bitstream;
  uint intprec = 32;
  uint kmin = intprec > maxprec ? intprec - maxprec : 0;

  // initialize data array to all zeros
  std::fill(data, data + size, (signed int)0);

  // input one bit plane at a time from MSB to LSB
  for (uint k = intprec, n = 0; k-- > kmin;) {
    // decode bit plane k
    for (uint i = 0;;) {
      // decode group of significant values
      for (; i < n; i++) {
        // read bit k of |data[i]|
        uint x = uint(stream.read()) << k;
        // NOTE: conditionals reordered to reduce branch mispredictions
        if (data[i])
          data[i] += data[i] < 0 ? -x : +x;
        else if (x) {
          // read sign bit also
          data[i] = stream.read() ? -x : +x;
        }
      }
      // have all groups been decoded?
      if (!count)
        break;
      // test next group
      if (stream.read()) {
        // group is significant; peel off and decode first subgroup
        n += count & 0xfu;
        count >>= 4;
      }
      else {
        // group is insignificant; continue with next bit plane
        break;
      }
    }
  }
exit:
  bitstream = stream;
}

#define ZFP_SZ_COMPRESS_FREQ 65536
#define ZFP_FACTOR_SAMPLE_NUM 8

int estimate_zfp_sz_compress_num(signed int * ordered_fp, int * block_exp, size_t num_sample_blocks, int block_size, double abs_eb, int zfp_factor){
	const int bs = 4;
    int sz_compress_num = 0;

	size_t dim0_offset = num_sample_blocks * block_size * block_size;
	size_t dim1_offset = num_sample_blocks * block_size;
	double coeff_eb = abs_eb * zfp_factor;
	int emin = INT_MIN;
	if (coeff_eb > 0) {
		frexp(coeff_eb, &emin);
		emin--;
	}
	double real_eb = ldexp(1.0, emin) / 64;
	int wlen = 32;
	int intbits = 2;

	size_t total_sample_count = 0;
	// ZFP statistics
	size_t encoding_length[9] = {0}; // encoding length for 9 groups
	size_t sign_bits[9] = {0};
	size_t group_bits[9] = {0};
	size_t skip_bits[9] = {0};
    int gp_element_num[9] = {1, 3, 6, 10, 12, 12, 10, 6, 4}; // number of data in a group

    // SZ statistics
    // size_t frequency[9][65536] = {0};
    size_t * frequency = (size_t *) malloc(ZFP_SZ_COMPRESS_FREQ*9*sizeof(size_t));
    memset(frequency, 0, ZFP_SZ_COMPRESS_FREQ*9*sizeof(size_t));
    int freq_radius = 32768;
    size_t unpred_count[9] = {0};

	float sample_coeff[64][2][2][2] = {0};
	signed int * ordered_fp_pos = ordered_fp;
	int * block_exp_pos = block_exp;
	for(size_t i=0; i<num_sample_blocks; i++){
		for(int ii=0; ii<2; ii++){
			for(int jj=0; jj<2; jj++){
				for(int kk=0; kk<2; kk++){
					// record coeff
					int e = *(block_exp_pos++);
					for(int num=0; num<64; num++){
						signed int fp = ordered_fp_pos[num];
				    	sample_coeff[num][ii][jj][kk] = ldexp((float)fp, intbits - wlen + e);
					}
					if(kk){
						int maxprec = MAX(0, e - emin + 8);
						uint intprec = 32;
						uint kmin = intprec > maxprec ? intprec - maxprec : 0;
						signed int * block_data_fp_pos = ordered_fp_pos;
						int last_max_ele = 32;
						for(int group=0; group<9; group++){
							int max_ele = 0;
							int cur_sign_bits = 0;
							for(int e=0; e<gp_element_num[group]; e++){
								uint fp = *block_data_fp_pos < 0 ? -uint(*block_data_fp_pos) : +uint(*block_data_fp_pos);
								int encoded_size = (int)fp_ufls(fp) - (int)kmin;
								if(encoded_size > 0) cur_sign_bits ++;
								max_ele = MAX(max_ele, encoded_size);
								block_data_fp_pos ++;
							}
							if(max_ele){
								encoding_length[group] += max_ele;
								group_bits[group] ++;
							}
							sign_bits[group] += cur_sign_bits;
							if(group > 0) skip_bits[group - 1] += last_max_ele - max_ele;
							last_max_ele = max_ele;
						}
						// estimate SZ efficiency
						int index = 0;
						for(int group=0; group<9; group++){
							for(int e=0; e<gp_element_num[group]; e++){
								float pred = sample_coeff[index][ii][jj][0];
								double diff = sample_coeff[index][ii][jj][1] - pred;
								double itvNum = fabs(diff)/real_eb + 1;
								if(itvNum < freq_radius){
									if(diff < 0) itvNum = -itvNum;
									frequency[group * ZFP_SZ_COMPRESS_FREQ + (int) (itvNum/2) + freq_radius] ++;
								}
								else{
									unpred_count[group] ++;
								}
								index ++;
							}
						}
						total_sample_count ++;						
					}
					ordered_fp_pos += 64;

				}
			}
		}
	}
	// sample done, estimate sz_compress_num
	// estimate zfp bit rate
	float zfp_br[9] = {0};
	float aver_zfp_br = 0;
	for(int i=0; i<9; i++){
		zfp_br[i] = encoding_length[i]*1.0/total_sample_count + group_bits[i]*1.0/total_sample_count/gp_element_num[i] + sign_bits[i]*1.0/total_sample_count/gp_element_num[i] + skip_bits[i]*1.0/total_sample_count/gp_element_num[i];
		aver_zfp_br += zfp_br[i] * gp_element_num[i];
	}
	aver_zfp_br /= 64;
	// estimate sz bit rate
	float sz_br[9] = {0};
	float aver_sz_br = 0;
	for(int i=0; i<9; i++){
		size_t encoded_num = total_sample_count*gp_element_num[i] - unpred_count[i];
		if(!encoded_num){
			sz_br[i] = 28;
			continue;
		}
		float br = 0;
		for(int e=0; e<ZFP_SZ_COMPRESS_FREQ; e++){
			if(frequency[i * ZFP_SZ_COMPRESS_FREQ + e]){
				float p = frequency[i * ZFP_SZ_COMPRESS_FREQ + e] * 1.0/encoded_num;
				br += -p * log2(p);
			}
		} 
		br += unpred_count[i] * 28.0 / (total_sample_count*gp_element_num[i]);
		sz_br[i] = br;
		aver_sz_br += br*gp_element_num[i];
	}
	aver_sz_br /= 64;
	free(frequency);

	sz_compress_num = 0;
	for(int i=0; i<9; i++){
		if(sz_br[i] < zfp_br[i]) sz_compress_num += gp_element_num[i];
		else break;
	}
	return sz_compress_num;
}

// compress on sampled data: block_size * block_size * (block_size * num_sample_blocks)
// given transformed coeffcients
void zfp_sample_compress(float * sample_data, signed int * sample_ordered_fp, int * block_exp, size_t num_sample_blocks, double abs_eb, int zfp_factor, int * zfp_sz_compress_num, double * mse, float * bit_rate){
	int block_size = 8;
	int sz_compress_num = estimate_zfp_sz_compress_num(sample_ordered_fp, block_exp, num_sample_blocks, block_size, abs_eb, zfp_factor);
	int sz_compress_group =0;
	int gp_element_num[9] = {1, 3, 6, 10, 12, 12, 10, 6, 4};
	int gp_element_offset[9] = {0, 1, 4, 10, 20, 32, 44, 54, 60};
	int64 zfp_modified_count;
	{
	    int64 count = 0x46acca631ull;
	    int tmp = sz_compress_num;
	    while(tmp > 0){
	    	tmp -= gp_element_num[sz_compress_group];
	    	count >>= 4;
	    	sz_compress_group ++;
	    }
	    zfp_modified_count = count;
	}
    int * combined_zfp_type = (int *) malloc(num_sample_blocks * block_size * block_size * block_size * sizeof(int));
    float * combined_zfp_unpred_data = (float *) malloc(num_sample_blocks * block_size * block_size * block_size * sizeof(float));
    unsigned char * combined_zfp_embedded = (unsigned char *) malloc(num_sample_blocks * block_size * block_size * block_size * sizeof(float));
    int * combined_zfp_sz_compress_coeff_type = (int *) malloc(num_sample_blocks * sz_compress_num * 8 * sizeof(int));
    float * combined_zfp_sz_compress_coeff_unpred_data = (float *) malloc(num_sample_blocks * sz_compress_num * 8 * sizeof(float));
    size_t combined_zfp_sz_compress_coeff_unpred_count = 0;
	size_t dim0_offset = num_sample_blocks * block_size * block_size;
	size_t dim1_offset = num_sample_blocks * block_size;
    ptrdiff_t zfp_offsets[8];
    for(int i=0; i<2; i++){
    	for(int j=0; j<2; j++){
    		for(int k=0; k<2; k++){
    			zfp_offsets[4*i+2*j+k] = i*4*dim0_offset + j*4*dim1_offset + k*4;
    		}
    	}
    }
    int zfp_order[8] = {6, 2, 0, 4, 5, 1, 3, 7};
    float * zfp_data_pos;
    float * zfp_last_sz_compressed_coeff = (float *) malloc(sz_compress_num*sizeof(float));
    memset(zfp_last_sz_compressed_coeff, 0, sz_compress_num * sizeof(float));

	int emin = INT_MIN;
	double realPrecision = abs_eb;
	double zfp_precision = realPrecision * zfp_factor;
	if (zfp_precision > 0) {
		frexp(zfp_precision, &emin);
		emin--;
	}
	double zfp_coeff_eb = ldexp(1.0, emin) / 64; // 4^3 = 64 is for the maximum loss in transform
	size_t zfp_sz_compress_coeff_index = 0;
	int * zfp_wrap_type = combined_zfp_type;
	float * zfp_wrap_unpred_data = combined_zfp_unpred_data;
	size_t zfp_wrap_unpred_count = 0;
	int zfp_coeff_intv_capacity = 65536;
	int zfp_coeff_intv_radius = 32768;
	int zfp_wrap_intv_capacity = 65536;
	int zfp_wrap_intv_radius = 32768;
	int wlen = 32;
	int intbits = 2;
	MemoryBitStream stream;
	stream.open(combined_zfp_embedded, num_sample_blocks*block_size*block_size*block_size*sizeof(float));
	double final_mse = 0;

	signed int block_data_fp[64];
	signed int tmp_ordered_fp[64];
	float block_data[64];
	float * sample_data_pos = sample_data;
	signed int * ordered_fp_pos = sample_ordered_fp;
	int * block_exp_pos = block_exp;
	for(size_t i=0; i<num_sample_blocks; i++){
		// compress by mZFP
		// each block includes 8 ZFP blocks
		for(int z=0; z<8; z++){
			signed int * cur_ordered_fp_pos = ordered_fp_pos + zfp_order[z]*64;
			memcpy(tmp_ordered_fp, ordered_fp_pos + zfp_order[z]*64, 64*sizeof(signed int));
			float * zfp_data_pos = sample_data_pos + zfp_offsets[zfp_order[z]];
			float max_ele = 0;
			for(int ii=0; ii<4; ii++){
				for(int jj=0; jj<4; jj++){
					for(int kk=0; kk<4; kk++){
						float cur_data = *(zfp_data_pos + ii*dim0_offset + jj*dim1_offset + kk);
						*(block_data + ii*16 + jj*4 + kk) = cur_data;
						max_ele = MAX(max_ele, fabs(cur_data));
					}
				}
			}
			// reorder exp because of zfp shift
			int e;
		    frexp(max_ele, &e);
		    *block_exp_pos = e;
			
			// coefficient analysis
			{
				memcpy(tmp_ordered_fp, cur_ordered_fp_pos, 64*sizeof(signed int));
				int e = *block_exp_pos;
				// compress sz_compress_coeff
				for(int num=0; num<sz_compress_num; num++){
					signed int fp = tmp_ordered_fp[num];
			    	float cur_data = ldexp((float)fp, intbits - wlen + e);
			    	// 
			    	float diff = cur_data - zfp_last_sz_compressed_coeff[num];
			    	double itvNum = fabs(diff) / zfp_coeff_eb + 1;
			    	if (itvNum < zfp_coeff_intv_capacity){
			    		if (diff < 0) itvNum = -itvNum;
			    		combined_zfp_sz_compress_coeff_type[zfp_sz_compress_coeff_index] = (int) (itvNum/2) + zfp_coeff_intv_radius;
						zfp_last_sz_compressed_coeff[num] = zfp_last_sz_compressed_coeff[num] + 2 * (combined_zfp_sz_compress_coeff_type[zfp_sz_compress_coeff_index] - zfp_coeff_intv_radius) * zfp_coeff_eb;
			    	}
			    	else{
			    		combined_zfp_sz_compress_coeff_type[zfp_sz_compress_coeff_index] = 0;
			    		zfp_last_sz_compressed_coeff[num] = cur_data;
			    		combined_zfp_sz_compress_coeff_unpred_data[combined_zfp_sz_compress_coeff_unpred_count ++] = cur_data;
			    	}
			    	zfp_sz_compress_coeff_index ++;
			    	// restore this coeff to fixed point
			    	cur_data = ldexp(zfp_last_sz_compressed_coeff[num], wlen - intbits - e);
			    	tmp_ordered_fp[num] = (long int) cur_data;
				}
				int maxprec = MAX(0, e - emin + 8);
				// embedded encoding
				uint m = 0;
				uint64 w = 0;
				for(int num=8; num>=sz_compress_group; num--){
					w = (w << 6) + fp_width(tmp_ordered_fp + gp_element_offset[num],  gp_element_num[num], m);
				}
				embedded_encoding(stream, tmp_ordered_fp + sz_compress_num, maxprec, zfp_modified_count, w);

				// restore block coeffcients
				uint intprec = 32;
				uint kmin = intprec > maxprec ? intprec - maxprec : 0;
				uint insignificant_bitplanes = kmin;
				for(int num=sz_compress_num; num<64; num++){
					signed int tmp_fp = (tmp_ordered_fp[num] > 0) ? tmp_ordered_fp[num] : - tmp_ordered_fp[num];
					// *******IMPORTANT REMEBER SIGN FOR ZFP***********
					signed int tmp = (tmp_fp >> insignificant_bitplanes) << insignificant_bitplanes;
					tmp = (tmp_ordered_fp[num] > 0) ? tmp : -tmp;
			    	tmp_ordered_fp[num] = tmp;
				}
				// reorder
				for(int ii=0; ii<4; ii++){
					for(int jj=0; jj<4; jj++){
						for(int kk=0; kk<4; kk++){
							int id = block_id[ii][jj][kk];
							*(block_data_fp + ii*16 + jj*4 + kk) = tmp_ordered_fp[id];
						}
					}
				}
			}
			// restore decompressed data
			// transform along x
			for(int jj=0; jj<4; jj++){
				for(int kk=0; kk<4; kk++){
					fp_inv_lift(block_data_fp + jj*4 + kk, 16);
				}
			}
			// transform along y
			for(int ii=0; ii<4; ii++){
				for(int kk=0; kk<4; kk++){
					fp_inv_lift(block_data_fp + ii*16 + kk, 4);
				}
			}
			// transform along z
			for(int ii=0; ii<4; ii++){
				for(int jj=0; jj<4; jj++){
					fp_inv_lift(block_data_fp + ii*16 + 4*jj, 1);
				}
			}
			// restore data && quantization
			{
				int e = *block_exp_pos;
				int wrap_index = 0;	
				for(int ii=0; ii<64; ii++){
			    	signed int fp = block_data_fp[ii];
			    	float pred = ldexp((float)fp, intbits - wlen + e);
			    	float diff = block_data[ii] - pred;
			    	double itvNum = fabs(diff) / realPrecision + 1;
			    	if (itvNum < zfp_wrap_intv_capacity){
			    		if (diff < 0) itvNum = -itvNum;
			    		zfp_wrap_type[wrap_index] = (int) (itvNum/2) + zfp_wrap_intv_radius;
			    		pred = pred + 2*(zfp_wrap_type[wrap_index] - zfp_wrap_intv_radius)*realPrecision;
			    		final_mse += (block_data[ii] - pred) * (block_data[ii] - pred);
			    	}
			    	else{
			    		zfp_wrap_type[wrap_index] = 0;
			    		zfp_wrap_unpred_data[zfp_wrap_unpred_count ++] = block_data[ii];
			    	}
			    	wrap_index ++;
			    }
			    zfp_wrap_type += 64;
			}
		    block_exp_pos ++;
		}
		ordered_fp_pos += 512;
		sample_data_pos += block_size;
	}
	stream.flush();
	unsigned char * comp_data = (unsigned char *) malloc(num_sample_blocks*block_size*block_size*block_size*sizeof(float));
	unsigned char * comp_data_pos = comp_data;
	SZ_Init(NULL);
	{
		*((float *) comp_data_pos) = zfp_factor;
		comp_data_pos += sizeof(float); 
		// block_exp
		size_t compressed_exp_size;
		unsigned char * compressed_exp = SZ_compress_args(SZ_INT32, block_exp, &compressed_exp_size, ABS, 0.5, 0, 0, 0, 0, 0, 0, num_sample_blocks*8);
		*((size_t *) comp_data_pos) = compressed_exp_size;
		comp_data_pos += sizeof(size_t);
		memcpy(comp_data_pos, compressed_exp, compressed_exp_size);
		comp_data_pos += compressed_exp_size;
		free(compressed_exp);
		// embedded encoding
		*((size_t *) comp_data_pos) = stream.size();
		comp_data_pos += sizeof(size_t);
		memcpy(comp_data_pos, combined_zfp_embedded, stream.size());
		comp_data_pos += stream.size();
		// zfp sz compressed coeff
		*((int *) comp_data_pos) = sz_compress_num;
		comp_data_pos += sizeof(int);
		if(sz_compress_num){
			unsigned char * tmp = comp_data_pos;
			// unpred data
			*((size_t *) comp_data_pos) = combined_zfp_sz_compress_coeff_unpred_count;
			comp_data_pos += sizeof(size_t);
			memcpy(comp_data_pos, combined_zfp_sz_compress_coeff_unpred_data, combined_zfp_sz_compress_coeff_unpred_count * sizeof(float));
			comp_data_pos += combined_zfp_sz_compress_coeff_unpred_count * sizeof(float);
			// huffman tree and type array
			int stateNum = zfp_coeff_intv_capacity;
			HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
			size_t nodeCount = 0;
			init(huffmanTree, combined_zfp_sz_compress_coeff_type, sz_compress_num*num_sample_blocks*8);
			size_t i = 0;
			for (i = 0; i < huffmanTree->stateNum; i++)
				if (huffmanTree->code[i]) nodeCount++; 
			nodeCount = nodeCount*2-1;
			unsigned char *treeBytes;
			unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(huffmanTree, nodeCount, &treeBytes);
			intToBytes_bigEndian(comp_data_pos, treeByteSize);
			comp_data_pos += sizeof(int);
			intToBytes_bigEndian(comp_data_pos, nodeCount);
			comp_data_pos += sizeof(int);
			memcpy(comp_data_pos, treeBytes, treeByteSize);
			comp_data_pos += treeByteSize;
			free(treeBytes);

			size_t typeArray_size = 0;
			encode(huffmanTree, combined_zfp_sz_compress_coeff_type, sz_compress_num*num_sample_blocks*8, comp_data_pos + sizeof(size_t), &typeArray_size);
			*((size_t *) comp_data_pos) = typeArray_size;
			comp_data_pos += sizeof(size_t) + typeArray_size;
			SZ_ReleaseHuffman(huffmanTree);
		}
		// wrapped
		{
			// put type array first for consistency
			int stateNum = zfp_wrap_intv_capacity;
			HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
			size_t nodeCount = 0;
			init(huffmanTree, combined_zfp_type, num_sample_blocks*block_size*block_size*block_size);
			size_t i = 0;
			for (i = 0; i < huffmanTree->stateNum; i++)
				if (huffmanTree->code[i]) nodeCount++; 
			nodeCount = nodeCount*2-1;
			unsigned char *treeBytes;
			unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(huffmanTree, nodeCount, &treeBytes);

			intToBytes_bigEndian(comp_data_pos, treeByteSize);
			comp_data_pos += sizeof(int);
			intToBytes_bigEndian(comp_data_pos, nodeCount);
			comp_data_pos += sizeof(int);
			memcpy(comp_data_pos, treeBytes, treeByteSize);
			comp_data_pos += treeByteSize;
			free(treeBytes);

			size_t typeArray_size = 0;
			encode(huffmanTree, combined_zfp_type, num_sample_blocks*block_size*block_size*block_size, comp_data_pos + sizeof(size_t), &typeArray_size);
			*((size_t *) comp_data_pos) = typeArray_size;
			comp_data_pos += sizeof(size_t) + typeArray_size;
			SZ_ReleaseHuffman(huffmanTree);
			// unpred data
			*((size_t *) comp_data_pos) = zfp_wrap_unpred_count;
			comp_data_pos += sizeof(size_t);
			memcpy(comp_data_pos, combined_zfp_unpred_data, zfp_wrap_unpred_count * sizeof(float));
			comp_data_pos += zfp_wrap_unpred_count * sizeof(float);
		}
	}
	unsigned char * compressed_comp_data = NULL;
	size_t totalEncodeSize = comp_data_pos - comp_data;
	size_t lossless_outsize = sz_lossless_compress(confparams_cpr->losslessCompressor, confparams_cpr->gzipMode, comp_data, totalEncodeSize, &compressed_comp_data);
	SZ_Finalize();
	free(compressed_comp_data);

	*zfp_sz_compress_num = sz_compress_num;
	*bit_rate = lossless_outsize * 32.0 / (num_sample_blocks * block_size * block_size * block_size * sizeof(float));
	*mse = final_mse/(num_sample_blocks * block_size * block_size * block_size);
	// printf("\nBit rate %.4g, mse %.4g\n\n", *bit_rate, *mse);
	free(comp_data);
	if(zfp_last_sz_compressed_coeff) free(zfp_last_sz_compressed_coeff);
	if(combined_zfp_embedded) free(combined_zfp_embedded);
	if(combined_zfp_type) free(combined_zfp_type);
	if(combined_zfp_unpred_data) free(combined_zfp_unpred_data);
	if(combined_zfp_sz_compress_coeff_type) free(combined_zfp_sz_compress_coeff_type);
	if(combined_zfp_sz_compress_coeff_unpred_data) free(combined_zfp_sz_compress_coeff_unpred_data);
}

unsigned char * zfp_compress(float * data, size_t r1, size_t r2, size_t r3, double abs_eb, int zfp_factor, int sz_compress_num, size_t * out_size, size_t * out_size_before_lossless, double * mse, float * bit_rate){
	int block_size = 8;
	size_t num_x = r1 / block_size;
	size_t num_y = r2 / block_size;
	size_t num_z = r3 / block_size;
	size_t num_blocks = num_x * num_y * num_z;
	size_t dim0_offset = r2*r3;
	size_t dim1_offset = r3;

	int sz_compress_group =0;
	int gp_element_num[9] = {1, 3, 6, 10, 12, 12, 10, 6, 4};
	int gp_element_offset[9] = {0, 1, 4, 10, 20, 32, 44, 54, 60};
	int64 zfp_modified_count;
	{
	    int64 count = 0x46acca631ull;
	    int tmp = sz_compress_num;
	    while(tmp > 0){
	    	tmp -= gp_element_num[sz_compress_group];
	    	count >>= 4;
	    	sz_compress_group ++;
	    }
	    zfp_modified_count = count;
	}
    int * combined_zfp_type = (int *) malloc(num_blocks * block_size * block_size * block_size * sizeof(int));
    float * combined_zfp_unpred_data = (float *) malloc(num_blocks * block_size * block_size * block_size * sizeof(float));
    unsigned char * combined_zfp_embedded = (unsigned char *) malloc(num_blocks * block_size * block_size * block_size * sizeof(float));
    int * combined_zfp_exp = (int *) malloc(num_blocks * 8 * sizeof(int));
    int * combined_zfp_sz_compress_coeff_type = (int *) malloc(num_blocks * sz_compress_num * 8 * sizeof(int));
    float * combined_zfp_sz_compress_coeff_unpred_data = (float *) malloc(num_blocks * sz_compress_num * 8 * sizeof(float));
    size_t combined_zfp_sz_compress_coeff_unpred_count = 0;
    ptrdiff_t zfp_offsets[8];
    for(int i=0; i<2; i++){
    	for(int j=0; j<2; j++){
    		for(int k=0; k<2; k++){
    			zfp_offsets[4*i+2*j+k] = i*4*dim0_offset + j*4*dim1_offset + k*4;
    		}
    	}
    }
    int zfp_order[8] = {6, 2, 0, 4, 5, 1, 3, 7};
    float * zfp_data_pos;
    float * zfp_last_sz_compressed_coeff = (float *) malloc(sz_compress_num*sizeof(float));
    memset(zfp_last_sz_compressed_coeff, 0, sz_compress_num * sizeof(float));
    int * block_exp_pos = combined_zfp_exp;

	int emin = INT_MIN;
	double realPrecision = abs_eb;
	double zfp_precision = realPrecision * zfp_factor;
	if (zfp_precision > 0) {
		frexp(zfp_precision, &emin);
		emin--;
	}
	double zfp_coeff_eb = ldexp(1.0, emin) / 64; // 4^3 = 64 is for the maximum loss in transform
	size_t zfp_sz_compress_coeff_index = 0;
	float block_data[64];
	signed int block_data_fp[64];
	signed int ordered_fp[64];
	int * zfp_wrap_type = combined_zfp_type;
	float * zfp_wrap_unpred_data = combined_zfp_unpred_data;
	size_t zfp_wrap_unpred_count = 0;
	int zfp_coeff_intv_capacity = 65536;
	int zfp_coeff_intv_radius = 32768;
	int zfp_wrap_intv_capacity = 65536;
	int zfp_wrap_intv_radius = 32768;
	int wlen = 32;
	int intbits = 2;
	MemoryBitStream stream;
	stream.open(combined_zfp_embedded, num_blocks*block_size*block_size*block_size*sizeof(float));
	double final_mse = 0;
    for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				float * data_pos = data + i*block_size * dim0_offset + j*block_size * dim1_offset + k*block_size;
				// compress by mZFP
				// each block includes 8 ZFP blocks
				for(int z=0; z<8; z++){
					float * zfp_data_pos = data_pos + zfp_offsets[zfp_order[z]];
					// encoding
					float max_ele = 0;
					for(int ii=0; ii<4; ii++){
						for(int jj=0; jj<4; jj++){
							for(int kk=0; kk<4; kk++){
								float cur_data = *(zfp_data_pos + ii*dim0_offset + jj*dim1_offset + kk);
								*(block_data + ii*16 + jj*4 + kk) = cur_data;
								max_ele = MAX(max_ele, fabs(cur_data));
							}
						}
					}
					// to fixed point
					{
						int e;
					    frexp(max_ele, &e);
					    *block_exp_pos = e;
						for(int ii=0; ii<64; ii++){
					    	float cur_data = ldexp(block_data[ii], wlen - intbits - e);
					    	block_data_fp[ii] = (long int) cur_data;
					    }
					}
					// transform
					// transform along z
					for(int ii=0; ii<4; ii++){
						for(int jj=0; jj<4; jj++){
							fp_fwd_lift(block_data_fp + ii*16 + 4*jj, 1);
						}
					}
					// transform along y
					for(int ii=0; ii<4; ii++){
						for(int kk=0; kk<4; kk++){
							fp_fwd_lift(block_data_fp + ii*16 + kk, 4);
						}
					}
					// transform along x
					for(int jj=0; jj<4; jj++){
						for(int kk=0; kk<4; kk++){
							fp_fwd_lift(block_data_fp + jj*4 + kk, 16);
						}
					}

					// coefficient analysis
					{
						// reorder
						for(int ii=0; ii<4; ii++){
							for(int jj=0; jj<4; jj++){
								for(int kk=0; kk<4; kk++){
									int id = block_id[ii][jj][kk];
									ordered_fp[id] = *(block_data_fp + ii*16 + jj*4 + kk);
								}
							}
						}
						int e = *block_exp_pos;
						// compress sz_compress_coeff
						for(int num=0; num<sz_compress_num; num++){
							signed int fp = ordered_fp[num];
					    	float cur_data = ldexp((float)fp, intbits - wlen + e);
					    	// 
					    	float diff = cur_data - zfp_last_sz_compressed_coeff[num];
					    	double itvNum = fabs(diff) / zfp_coeff_eb + 1;
					    	if (itvNum < zfp_coeff_intv_capacity){
					    		if (diff < 0) itvNum = -itvNum;
					    		combined_zfp_sz_compress_coeff_type[zfp_sz_compress_coeff_index] = (int) (itvNum/2) + zfp_coeff_intv_radius;
								zfp_last_sz_compressed_coeff[num] = zfp_last_sz_compressed_coeff[num] + 2 * (combined_zfp_sz_compress_coeff_type[zfp_sz_compress_coeff_index] - zfp_coeff_intv_radius) * zfp_coeff_eb;
					    	}
					    	else{
					    		combined_zfp_sz_compress_coeff_type[zfp_sz_compress_coeff_index] = 0;
					    		zfp_last_sz_compressed_coeff[num] = cur_data;
					    		combined_zfp_sz_compress_coeff_unpred_data[combined_zfp_sz_compress_coeff_unpred_count ++] = cur_data;
					    	}
					    	zfp_sz_compress_coeff_index ++;
					    	// restore this coeff to fixed point
					    	cur_data = ldexp(zfp_last_sz_compressed_coeff[num], wlen - intbits - e);
					    	ordered_fp[num] = (long int) cur_data;
						}
						int maxprec = MAX(0, e - emin + 8);
						// embedded encoding
						uint m = 0;
						uint64 w = 0;
						for(int num=8; num>=sz_compress_group; num--){
							w = (w << 6) + fp_width(ordered_fp + gp_element_offset[num],  gp_element_num[num], m);
						}
						embedded_encoding(stream, ordered_fp + sz_compress_num, maxprec, zfp_modified_count, w);

						// restore block coeffcients
						uint intprec = 32;
						uint kmin = intprec > maxprec ? intprec - maxprec : 0;
						uint insignificant_bitplanes = kmin;
						for(int num=sz_compress_num; num<64; num++){
							// ordered_fp[num] = (ordered_fp[num] >> insignificant_bitplanes) << insignificant_bitplanes;
							signed int tmp_fp = (ordered_fp[num] > 0) ? ordered_fp[num] : - ordered_fp[num];
							// *******IMPORTANT REMEBER SIGN FOR ZFP***********
							signed int tmp = (tmp_fp >> insignificant_bitplanes) << insignificant_bitplanes;
							tmp = (ordered_fp[num] > 0) ? tmp : -tmp;
					    	ordered_fp[num] = tmp;
						}
						// reorder
						for(int ii=0; ii<4; ii++){
							for(int jj=0; jj<4; jj++){
								for(int kk=0; kk<4; kk++){
									int id = block_id[ii][jj][kk];
									*(block_data_fp + ii*16 + jj*4 + kk) = ordered_fp[id];
								}
							}
						}
					}
					// restore decompressed data
					// transform along x
					for(int jj=0; jj<4; jj++){
						for(int kk=0; kk<4; kk++){
							fp_inv_lift(block_data_fp + jj*4 + kk, 16);
						}
					}
					// transform along y
					for(int ii=0; ii<4; ii++){
						for(int kk=0; kk<4; kk++){
							fp_inv_lift(block_data_fp + ii*16 + kk, 4);
						}
					}
					// transform along z
					for(int ii=0; ii<4; ii++){
						for(int jj=0; jj<4; jj++){
							fp_inv_lift(block_data_fp + ii*16 + 4*jj, 1);
						}
					}
					// restore data && quantization
					{
						int e = *block_exp_pos;
						int wrap_index = 0;	
						for(int ii=0; ii<64; ii++){
					    	signed int fp = block_data_fp[ii];
					    	float pred = ldexp((float)fp, intbits - wlen + e);
					    	float diff = block_data[ii] - pred;
					    	double itvNum = fabs(diff) / realPrecision + 1;
					    	if (itvNum < zfp_wrap_intv_capacity){
					    		if (diff < 0) itvNum = -itvNum;
					    		zfp_wrap_type[wrap_index] = (int) (itvNum/2) + zfp_wrap_intv_radius;
					    		pred = pred + 2*(zfp_wrap_type[wrap_index] - zfp_wrap_intv_radius)*realPrecision;
					    		final_mse += (block_data[ii] - pred) * (block_data[ii] - pred);
					    	}
					    	else{
					    		zfp_wrap_type[wrap_index] = 0;
					    		zfp_wrap_unpred_data[zfp_wrap_unpred_count ++] = block_data[ii];
					    	}
					    	wrap_index ++;
					    }
					    zfp_wrap_type += 64;
					}
				    block_exp_pos ++;
				}
			}
		}
	}
	stream.flush();
	unsigned char * comp_data = (unsigned char *) malloc(num_blocks*block_size*block_size*block_size*sizeof(float));
	unsigned char * comp_data_pos = comp_data;
	SZ_Init(NULL);
	{
		*((double *) comp_data_pos) = realPrecision;
		comp_data_pos += sizeof(double); 
		*((int *) comp_data_pos) = zfp_factor;
		comp_data_pos += sizeof(float); 
		// block_exp
		size_t compressed_exp_size;
		unsigned char * compressed_exp = SZ_compress_args(SZ_INT32, combined_zfp_exp, &compressed_exp_size, ABS, 0.5, 0, 0, 0, 0, 0, 0, num_blocks*8);
		*((size_t *) comp_data_pos) = compressed_exp_size;
		comp_data_pos += sizeof(size_t);
		memcpy(comp_data_pos, compressed_exp, compressed_exp_size);
		comp_data_pos += compressed_exp_size;
		free(compressed_exp);
		// embedded encoding
		*((size_t *) comp_data_pos) = stream.size();
		comp_data_pos += sizeof(size_t);
		memcpy(comp_data_pos, combined_zfp_embedded, stream.size());
		comp_data_pos += stream.size();
		// zfp sz compressed coeff
		*((int *) comp_data_pos) = sz_compress_num;
		comp_data_pos += sizeof(int);
		if(sz_compress_num){
			unsigned char * tmp = comp_data_pos;
			// unpred data
			*((size_t *) comp_data_pos) = combined_zfp_sz_compress_coeff_unpred_count;
			comp_data_pos += sizeof(size_t);
			memcpy(comp_data_pos, combined_zfp_sz_compress_coeff_unpred_data, combined_zfp_sz_compress_coeff_unpred_count * sizeof(float));
			comp_data_pos += combined_zfp_sz_compress_coeff_unpred_count * sizeof(float);
			// huffman tree and type array
			int stateNum = zfp_coeff_intv_capacity;
			HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
			size_t nodeCount = 0;
			init(huffmanTree, combined_zfp_sz_compress_coeff_type, sz_compress_num*num_blocks*8);
			size_t i = 0;
			for (i = 0; i < huffmanTree->stateNum; i++)
				if (huffmanTree->code[i]) nodeCount++; 
			nodeCount = nodeCount*2-1;
			unsigned char *treeBytes;
			unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(huffmanTree, nodeCount, &treeBytes);
			intToBytes_bigEndian(comp_data_pos, treeByteSize);
			comp_data_pos += sizeof(int);
			intToBytes_bigEndian(comp_data_pos, nodeCount);
			comp_data_pos += sizeof(int);
			memcpy(comp_data_pos, treeBytes, treeByteSize);
			comp_data_pos += treeByteSize;
			free(treeBytes);

			size_t typeArray_size = 0;
			encode(huffmanTree, combined_zfp_sz_compress_coeff_type, sz_compress_num*num_blocks*8, comp_data_pos + sizeof(size_t), &typeArray_size);
			*((size_t *) comp_data_pos) = typeArray_size;
			comp_data_pos += sizeof(size_t) + typeArray_size;
			SZ_ReleaseHuffman(huffmanTree);
		}
		// wrapped
		{
			// put type array first for consistency
			int stateNum = zfp_wrap_intv_capacity;
			HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
			size_t nodeCount = 0;
			init(huffmanTree, combined_zfp_type, num_blocks*block_size*block_size*block_size);
			size_t i = 0;
			for (i = 0; i < huffmanTree->stateNum; i++)
				if (huffmanTree->code[i]) nodeCount++; 
			nodeCount = nodeCount*2-1;
			unsigned char *treeBytes;
			unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(huffmanTree, nodeCount, &treeBytes);

			intToBytes_bigEndian(comp_data_pos, treeByteSize);
			comp_data_pos += sizeof(int);
			intToBytes_bigEndian(comp_data_pos, nodeCount);
			comp_data_pos += sizeof(int);
			memcpy(comp_data_pos, treeBytes, treeByteSize);
			comp_data_pos += treeByteSize;
			free(treeBytes);

			size_t typeArray_size = 0;
			encode(huffmanTree, combined_zfp_type, num_blocks*block_size*block_size*block_size, comp_data_pos + sizeof(size_t), &typeArray_size);
			*((size_t *) comp_data_pos) = typeArray_size;
			comp_data_pos += sizeof(size_t) + typeArray_size;
			SZ_ReleaseHuffman(huffmanTree);
			// unpred data
			*((size_t *) comp_data_pos) = zfp_wrap_unpred_count;
			comp_data_pos += sizeof(size_t);
			memcpy(comp_data_pos, combined_zfp_unpred_data, zfp_wrap_unpred_count * sizeof(float));
			comp_data_pos += zfp_wrap_unpred_count * sizeof(float);
		}
	}
	unsigned char * compressed_comp_data = NULL;
	size_t totalEncodeSize = comp_data_pos - comp_data;
	size_t lossless_outsize = sz_lossless_compress(confparams_cpr->losslessCompressor, confparams_cpr->gzipMode, comp_data, totalEncodeSize, &compressed_comp_data);
	SZ_Finalize();
	*out_size = lossless_outsize;
	*out_size_before_lossless = totalEncodeSize;
	*bit_rate = lossless_outsize * 32.0 / (num_blocks * block_size * block_size * block_size * sizeof(float));
	*mse = final_mse/(num_blocks * block_size * block_size * block_size);
	free(comp_data);
	if(zfp_last_sz_compressed_coeff) free(zfp_last_sz_compressed_coeff);
	if(combined_zfp_exp) free(combined_zfp_exp);
	if(combined_zfp_embedded) free(combined_zfp_embedded);
	if(combined_zfp_type) free(combined_zfp_type);
	if(combined_zfp_unpred_data) free(combined_zfp_unpred_data);
	if(combined_zfp_sz_compress_coeff_type) free(combined_zfp_sz_compress_coeff_type);
	if(combined_zfp_sz_compress_coeff_unpred_data) free(combined_zfp_sz_compress_coeff_unpred_data);
	return compressed_comp_data;
}

unsigned char * zfp_compress_block_size_4(float * data, size_t r1, size_t r2, size_t r3, double abs_eb, int zfp_factor, int sz_compress_num, size_t * out_size, size_t * out_size_before_lossless, double * mse, float * bit_rate){
	int block_size = 4;
	size_t num_x = r1 / block_size;
	size_t num_y = r2 / block_size;
	size_t num_z = r3 / block_size;
	size_t num_blocks = num_x * num_y * num_z;
	size_t dim0_offset = r2*r3;
	size_t dim1_offset = r3;

	int sz_compress_group =0;
	int gp_element_num[9] = {1, 3, 6, 10, 12, 12, 10, 6, 4};
	int gp_element_offset[9] = {0, 1, 4, 10, 20, 32, 44, 54, 60};
	int64 zfp_modified_count;
	{
	    int64 count = 0x46acca631ull;
	    int tmp = sz_compress_num;
	    while(tmp > 0){
	    	tmp -= gp_element_num[sz_compress_group];
	    	count >>= 4;
	    	sz_compress_group ++;
	    }
	    zfp_modified_count = count;
	}
    int * combined_zfp_type = (int *) malloc(num_blocks * block_size * block_size * block_size * sizeof(int));
    float * combined_zfp_unpred_data = (float *) malloc(num_blocks * block_size * block_size * block_size * sizeof(float));
    unsigned char * combined_zfp_embedded = (unsigned char *) malloc(num_blocks * block_size * block_size * block_size * sizeof(float));
    int * combined_zfp_exp = (int *) malloc(num_blocks * sizeof(int));
    int * combined_zfp_sz_compress_coeff_type = (int *) malloc(num_blocks * sz_compress_num * sizeof(int));
    float * combined_zfp_sz_compress_coeff_unpred_data = (float *) malloc(num_blocks * sz_compress_num * sizeof(float));
    size_t combined_zfp_sz_compress_coeff_unpred_count = 0;

    float * zfp_last_sz_compressed_coeff = (float *) malloc(sz_compress_num*sizeof(float));
    memset(zfp_last_sz_compressed_coeff, 0, sz_compress_num * sizeof(float));
    int * block_exp_pos = combined_zfp_exp;

	int emin = INT_MIN;
	double realPrecision = abs_eb;
	double zfp_precision = realPrecision * zfp_factor;
	if (zfp_precision > 0) {
		frexp(zfp_precision, &emin);
		emin--;
	}
	double zfp_coeff_eb = ldexp(1.0, emin) / 64; // 4^3 = 64 is for the maximum loss in transform
	size_t zfp_sz_compress_coeff_index = 0;
	float block_data[64];
	signed int block_data_fp[64];
	signed int ordered_fp[64];
	int * zfp_wrap_type = combined_zfp_type;
	float * zfp_wrap_unpred_data = combined_zfp_unpred_data;
	size_t zfp_wrap_unpred_count = 0;
	int zfp_coeff_intv_capacity = 65536;
	int zfp_coeff_intv_radius = 32768;
	int zfp_wrap_intv_capacity = 65536;
	int zfp_wrap_intv_radius = 32768;
	int wlen = 32;
	int intbits = 2;
	MemoryBitStream stream;
	stream.open(combined_zfp_embedded, num_blocks*block_size*block_size*block_size*sizeof(float));
	double final_mse = 0;
    for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				float * data_pos = data + i*block_size * dim0_offset + j*block_size * dim1_offset + k*block_size;
				// compress by mZFP
				float * zfp_data_pos = data_pos;
				// encoding
				float max_ele = 0;
				for(int ii=0; ii<4; ii++){
					for(int jj=0; jj<4; jj++){
						for(int kk=0; kk<4; kk++){
							float cur_data = *(zfp_data_pos + ii*dim0_offset + jj*dim1_offset + kk);
							*(block_data + ii*16 + jj*4 + kk) = cur_data;
							max_ele = MAX(max_ele, fabs(cur_data));
						}
					}
				}
				// to fixed point
				{
					int e;
				    frexp(max_ele, &e);
				    *block_exp_pos = e;
					for(int ii=0; ii<64; ii++){
				    	float cur_data = ldexp(block_data[ii], wlen - intbits - e);
				    	block_data_fp[ii] = (long int) cur_data;
				    }
				}
				// transform
				// transform along z
				for(int ii=0; ii<4; ii++){
					for(int jj=0; jj<4; jj++){
						fp_fwd_lift(block_data_fp + ii*16 + 4*jj, 1);
					}
				}
				// transform along y
				for(int ii=0; ii<4; ii++){
					for(int kk=0; kk<4; kk++){
						fp_fwd_lift(block_data_fp + ii*16 + kk, 4);
					}
				}
				// transform along x
				for(int jj=0; jj<4; jj++){
					for(int kk=0; kk<4; kk++){
						fp_fwd_lift(block_data_fp + jj*4 + kk, 16);
					}
				}

				// coefficient analysis
				{
					// reorder
					for(int ii=0; ii<4; ii++){
						for(int jj=0; jj<4; jj++){
							for(int kk=0; kk<4; kk++){
								int id = block_id[ii][jj][kk];
								ordered_fp[id] = *(block_data_fp + ii*16 + jj*4 + kk);
							}
						}
					}
					int e = *block_exp_pos;
					// compress sz_compress_coeff
					for(int num=0; num<sz_compress_num; num++){
						signed int fp = ordered_fp[num];
				    	float cur_data = ldexp((float)fp, intbits - wlen + e);
				    	// 
				    	float diff = cur_data - zfp_last_sz_compressed_coeff[num];
				    	double itvNum = fabs(diff) / zfp_coeff_eb + 1;
				    	if (itvNum < zfp_coeff_intv_capacity){
				    		if (diff < 0) itvNum = -itvNum;
				    		combined_zfp_sz_compress_coeff_type[zfp_sz_compress_coeff_index] = (int) (itvNum/2) + zfp_coeff_intv_radius;
							zfp_last_sz_compressed_coeff[num] = zfp_last_sz_compressed_coeff[num] + 2 * (combined_zfp_sz_compress_coeff_type[zfp_sz_compress_coeff_index] - zfp_coeff_intv_radius) * zfp_coeff_eb;
				    	}
				    	else{
				    		combined_zfp_sz_compress_coeff_type[zfp_sz_compress_coeff_index] = 0;
				    		zfp_last_sz_compressed_coeff[num] = cur_data;
				    		combined_zfp_sz_compress_coeff_unpred_data[combined_zfp_sz_compress_coeff_unpred_count ++] = cur_data;
				    	}
				    	zfp_sz_compress_coeff_index ++;
				    	// restore this coeff to fixed point
				    	cur_data = ldexp(zfp_last_sz_compressed_coeff[num], wlen - intbits - e);
				    	ordered_fp[num] = (long int) cur_data;
					}
					int maxprec = MAX(0, e - emin + 8);
					// embedded encoding
					uint m = 0;
					uint64 w = 0;
					for(int num=8; num>=sz_compress_group; num--){
						w = (w << 6) + fp_width(ordered_fp + gp_element_offset[num],  gp_element_num[num], m);
					}
					embedded_encoding(stream, ordered_fp + sz_compress_num, maxprec, zfp_modified_count, w);

					// restore block coeffcients
					uint intprec = 32;
					uint kmin = intprec > maxprec ? intprec - maxprec : 0;
					uint insignificant_bitplanes = kmin;
					for(int num=sz_compress_num; num<64; num++){
						// ordered_fp[num] = (ordered_fp[num] >> insignificant_bitplanes) << insignificant_bitplanes;
						signed int tmp_fp = (ordered_fp[num] > 0) ? ordered_fp[num] : - ordered_fp[num];
						// *******IMPORTANT REMEBER SIGN FOR ZFP***********
						signed int tmp = (tmp_fp >> insignificant_bitplanes) << insignificant_bitplanes;
						tmp = (ordered_fp[num] > 0) ? tmp : -tmp;
				    	ordered_fp[num] = tmp;
					}
					// reorder
					for(int ii=0; ii<4; ii++){
						for(int jj=0; jj<4; jj++){
							for(int kk=0; kk<4; kk++){
								int id = block_id[ii][jj][kk];
								*(block_data_fp + ii*16 + jj*4 + kk) = ordered_fp[id];
							}
						}
					}
				}
				// restore decompressed data
				// transform along x
				for(int jj=0; jj<4; jj++){
					for(int kk=0; kk<4; kk++){
						fp_inv_lift(block_data_fp + jj*4 + kk, 16);
					}
				}
				// transform along y
				for(int ii=0; ii<4; ii++){
					for(int kk=0; kk<4; kk++){
						fp_inv_lift(block_data_fp + ii*16 + kk, 4);
					}
				}
				// transform along z
				for(int ii=0; ii<4; ii++){
					for(int jj=0; jj<4; jj++){
						fp_inv_lift(block_data_fp + ii*16 + 4*jj, 1);
					}
				}
				// restore data && quantization
				{
					int e = *block_exp_pos;
					int wrap_index = 0;	
					for(int ii=0; ii<64; ii++){
				    	signed int fp = block_data_fp[ii];
				    	float pred = ldexp((float)fp, intbits - wlen + e);
				    	float diff = block_data[ii] - pred;
				    	double itvNum = fabs(diff) / realPrecision + 1;
				    	if (itvNum < zfp_wrap_intv_capacity){
				    		if (diff < 0) itvNum = -itvNum;
				    		zfp_wrap_type[wrap_index] = (int) (itvNum/2) + zfp_wrap_intv_radius;
				    		pred = pred + 2*(zfp_wrap_type[wrap_index] - zfp_wrap_intv_radius)*realPrecision;
				    		final_mse += (block_data[ii] - pred) * (block_data[ii] - pred);
				    	}
				    	else{
				    		zfp_wrap_type[wrap_index] = 0;
				    		zfp_wrap_unpred_data[zfp_wrap_unpred_count ++] = block_data[ii];
				    	}
				    	wrap_index ++;
				    }
				    zfp_wrap_type += 64;
				}
			    block_exp_pos ++;
			}
		}
	}
	stream.flush();
	unsigned char * comp_data = (unsigned char *) malloc(num_blocks*block_size*block_size*block_size*sizeof(float));
	unsigned char * comp_data_pos = comp_data;
	SZ_Init(NULL);
	{
		*((double *) comp_data_pos) = realPrecision;
		comp_data_pos += sizeof(double); 
		*((int *) comp_data_pos) = zfp_factor;
		comp_data_pos += sizeof(float); 
		// block_exp
		size_t compressed_exp_size;
		unsigned char * compressed_exp = SZ_compress_args(SZ_INT32, combined_zfp_exp, &compressed_exp_size, ABS, 0.5, 0, 0, 0, 0, 0, 0, num_blocks);
		*((size_t *) comp_data_pos) = compressed_exp_size;
		comp_data_pos += sizeof(size_t);
		memcpy(comp_data_pos, compressed_exp, compressed_exp_size);
		comp_data_pos += compressed_exp_size;
		free(compressed_exp);
		// embedded encoding
		*((size_t *) comp_data_pos) = stream.size();
		comp_data_pos += sizeof(size_t);
		memcpy(comp_data_pos, combined_zfp_embedded, stream.size());
		comp_data_pos += stream.size();
		// zfp sz compressed coeff
		*((int *) comp_data_pos) = sz_compress_num;
		comp_data_pos += sizeof(int);
		if(sz_compress_num){
			unsigned char * tmp = comp_data_pos;
			// unpred data
			*((size_t *) comp_data_pos) = combined_zfp_sz_compress_coeff_unpred_count;
			comp_data_pos += sizeof(size_t);
			memcpy(comp_data_pos, combined_zfp_sz_compress_coeff_unpred_data, combined_zfp_sz_compress_coeff_unpred_count * sizeof(float));
			comp_data_pos += combined_zfp_sz_compress_coeff_unpred_count * sizeof(float);
			// huffman tree and type array
			int stateNum = zfp_coeff_intv_capacity;
			HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
			size_t nodeCount = 0;
			init(huffmanTree, combined_zfp_sz_compress_coeff_type, sz_compress_num*num_blocks);
			size_t i = 0;
			for (i = 0; i < huffmanTree->stateNum; i++)
				if (huffmanTree->code[i]) nodeCount++; 
			nodeCount = nodeCount*2-1;
			unsigned char *treeBytes;
			unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(huffmanTree, nodeCount, &treeBytes);
			intToBytes_bigEndian(comp_data_pos, treeByteSize);
			comp_data_pos += sizeof(int);
			intToBytes_bigEndian(comp_data_pos, nodeCount);
			comp_data_pos += sizeof(int);
			memcpy(comp_data_pos, treeBytes, treeByteSize);
			comp_data_pos += treeByteSize;
			free(treeBytes);

			size_t typeArray_size = 0;
			encode(huffmanTree, combined_zfp_sz_compress_coeff_type, sz_compress_num*num_blocks, comp_data_pos + sizeof(size_t), &typeArray_size);
			*((size_t *) comp_data_pos) = typeArray_size;
			comp_data_pos += sizeof(size_t) + typeArray_size;
			SZ_ReleaseHuffman(huffmanTree);
		}
		// wrapped
		{
			// put type array first for consistency
			int stateNum = zfp_wrap_intv_capacity;
			HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
			size_t nodeCount = 0;
			init(huffmanTree, combined_zfp_type, num_blocks*block_size*block_size*block_size);
			size_t i = 0;
			for (i = 0; i < huffmanTree->stateNum; i++)
				if (huffmanTree->code[i]) nodeCount++; 
			nodeCount = nodeCount*2-1;
			unsigned char *treeBytes;
			unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(huffmanTree, nodeCount, &treeBytes);

			intToBytes_bigEndian(comp_data_pos, treeByteSize);
			comp_data_pos += sizeof(int);
			intToBytes_bigEndian(comp_data_pos, nodeCount);
			comp_data_pos += sizeof(int);
			memcpy(comp_data_pos, treeBytes, treeByteSize);
			comp_data_pos += treeByteSize;
			free(treeBytes);

			size_t typeArray_size = 0;
			encode(huffmanTree, combined_zfp_type, num_blocks*block_size*block_size*block_size, comp_data_pos + sizeof(size_t), &typeArray_size);
			*((size_t *) comp_data_pos) = typeArray_size;
			comp_data_pos += sizeof(size_t) + typeArray_size;
			SZ_ReleaseHuffman(huffmanTree);
			// unpred data
			*((size_t *) comp_data_pos) = zfp_wrap_unpred_count;
			comp_data_pos += sizeof(size_t);
			memcpy(comp_data_pos, combined_zfp_unpred_data, zfp_wrap_unpred_count * sizeof(float));
			comp_data_pos += zfp_wrap_unpred_count * sizeof(float);
		}
	}
	unsigned char * compressed_comp_data = NULL;
	size_t totalEncodeSize = comp_data_pos - comp_data;
	size_t lossless_outsize = sz_lossless_compress(confparams_cpr->losslessCompressor, confparams_cpr->gzipMode, comp_data, totalEncodeSize, &compressed_comp_data);
	SZ_Finalize();
	*out_size = lossless_outsize;
	*out_size_before_lossless = totalEncodeSize;
	*bit_rate = lossless_outsize * 32.0 / (num_blocks * block_size * block_size * block_size * sizeof(float));
	*mse = final_mse/(num_blocks * block_size * block_size * block_size);
	free(comp_data);
	if(zfp_last_sz_compressed_coeff) free(zfp_last_sz_compressed_coeff);
	if(combined_zfp_exp) free(combined_zfp_exp);
	if(combined_zfp_embedded) free(combined_zfp_embedded);
	if(combined_zfp_type) free(combined_zfp_type);
	if(combined_zfp_unpred_data) free(combined_zfp_unpred_data);
	if(combined_zfp_sz_compress_coeff_type) free(combined_zfp_sz_compress_coeff_type);
	if(combined_zfp_sz_compress_coeff_unpred_data) free(combined_zfp_sz_compress_coeff_unpred_data);
	return compressed_comp_data;
}

void sz_compute_sample_reg_coefficient(float * sample_data, size_t num_sample_blocks, int block_size, int use_mean, float mean, float * reg_params, float * reg_pred_err, float * sz_pred_err){

	size_t dim0_offset = num_sample_blocks * block_size * block_size;
	size_t dim1_offset = num_sample_blocks * block_size;
	int params_offset_b = num_sample_blocks;
	int params_offset_c = 2*num_sample_blocks;
	int params_offset_d = 3*num_sample_blocks;

	float * reg_params_pos = reg_params;
	size_t pred_err_index = 0;
	float * sample_data_pos = sample_data;
	for(size_t i=0; i<num_sample_blocks; i++){
		/*Calculate regression coefficients*/
		{
			float * cur_data_pos = sample_data_pos;
			float fx = 0.0;
			float fy = 0.0;
			float fz = 0.0;
			float f = 0;
			float sum_x, sum_y; 
			float curData;
			for(size_t i=0; i<block_size; i++){
				sum_x = 0;
				for(size_t j=0; j<block_size; j++){
					sum_y = 0;
					for(size_t k=0; k<block_size; k++){
						curData = *cur_data_pos;
						sum_y += curData;
						fz += curData * k;
						cur_data_pos ++;
					}
					fy += sum_y * j;
					sum_x += sum_y;
					cur_data_pos += dim1_offset - block_size;
				}
				fx += sum_x * i;
				f += sum_x;
				cur_data_pos += dim0_offset - block_size * dim1_offset;
			}
			float coeff = 1.0 / (block_size * block_size * block_size);
			reg_params_pos[0] = (2 * fx / (block_size - 1) - f) * 6 * coeff / (block_size + 1);
			reg_params_pos[params_offset_b] = (2 * fy / (block_size - 1) - f) * 6 * coeff / (block_size + 1);
			reg_params_pos[params_offset_c] = (2 * fz / (block_size - 1) - f) * 6 * coeff / (block_size + 1);
			reg_params_pos[params_offset_d] = f * coeff - ((block_size - 1) * reg_params_pos[0] / 2 + (block_size - 1) * reg_params_pos[params_offset_b] / 2 + (block_size - 1) * reg_params_pos[params_offset_c] / 2);
		}
		// decide
		/*sampling and decide which predictor*/
		// no noise here, will add when decide
		if(use_mean){
			// sample point [1, 1, 1] [1, 1, 4] [1, 4, 1] [1, 4, 4] [4, 1, 1] [4, 1, 4] [4, 4, 1] [4, 4, 4]
			float * data_pos = sample_data_pos;
			float * cur_data_pos;
			float curData;
			float pred_reg, pred_sz;
			float err_sz = 0.0, err_reg = 0.0;
			int bmi = 0;
			for(int i=1; i<block_size; i++){
				cur_data_pos = data_pos + i*dim0_offset + i*dim1_offset + i;
				curData = *cur_data_pos;
				pred_sz = cur_data_pos[-1] + cur_data_pos[-dim1_offset]+ cur_data_pos[-dim0_offset] - cur_data_pos[-dim1_offset - 1] - cur_data_pos[-dim0_offset - 1] - cur_data_pos[-dim0_offset - dim1_offset] + cur_data_pos[-dim0_offset - dim1_offset - 1];
				pred_reg = reg_params_pos[0] * i + reg_params_pos[params_offset_b] * i + reg_params_pos[params_offset_c] * i + reg_params_pos[params_offset_d];							
				err_sz += MIN(fabs(pred_sz - curData), fabs(mean - curData));
				err_reg += fabs(pred_reg - curData);

				bmi = block_size - i;
				cur_data_pos = data_pos + i*dim0_offset + i*dim1_offset + bmi;
				curData = *cur_data_pos;
				pred_sz = cur_data_pos[-1] + cur_data_pos[-dim1_offset]+ cur_data_pos[-dim0_offset] - cur_data_pos[-dim1_offset - 1] - cur_data_pos[-dim0_offset - 1] - cur_data_pos[-dim0_offset - dim1_offset] + cur_data_pos[-dim0_offset - dim1_offset - 1];
				pred_reg = reg_params_pos[0] * i + reg_params_pos[params_offset_b] * i + reg_params_pos[params_offset_c] * bmi + reg_params_pos[params_offset_d];							
				err_sz += MIN(fabs(pred_sz - curData), fabs(mean - curData));
				err_reg += fabs(pred_reg - curData);								

				cur_data_pos = data_pos + i*dim0_offset + bmi*dim1_offset + i;
				curData = *cur_data_pos;
				pred_sz = cur_data_pos[-1] + cur_data_pos[-dim1_offset]+ cur_data_pos[-dim0_offset] - cur_data_pos[-dim1_offset - 1] - cur_data_pos[-dim0_offset - 1] - cur_data_pos[-dim0_offset - dim1_offset] + cur_data_pos[-dim0_offset - dim1_offset - 1];
				pred_reg = reg_params_pos[0] * i + reg_params_pos[params_offset_b] * bmi + reg_params_pos[params_offset_c] * i + reg_params_pos[params_offset_d];							
				err_sz += MIN(fabs(pred_sz - curData), fabs(mean - curData));
				err_reg += fabs(pred_reg - curData);								

				cur_data_pos = data_pos + i*dim0_offset + bmi*dim1_offset + bmi;
				curData = *cur_data_pos;
				pred_sz = cur_data_pos[-1] + cur_data_pos[-dim1_offset]+ cur_data_pos[-dim0_offset] - cur_data_pos[-dim1_offset - 1] - cur_data_pos[-dim0_offset - 1] - cur_data_pos[-dim0_offset - dim1_offset] + cur_data_pos[-dim0_offset - dim1_offset - 1];
				pred_reg = reg_params_pos[0] * i + reg_params_pos[params_offset_b] * bmi + reg_params_pos[params_offset_c] * bmi + reg_params_pos[params_offset_d];							
				err_sz += MIN(fabs(pred_sz - curData), fabs(mean - curData));
				err_reg += fabs(pred_reg - curData);								

			}
			reg_pred_err[pred_err_index] = err_reg;
			sz_pred_err[pred_err_index] = err_sz;
			pred_err_index ++;
		}
		else{
			float * data_pos = sample_data_pos;
			float * cur_data_pos;
			float curData;
			float pred_reg, pred_sz;
			float err_sz = 0.0, err_reg = 0.0;
			int bmi = 0;							
			for(int i=1; i<block_size; i++){
				cur_data_pos = data_pos + i*dim0_offset + i*dim1_offset + i;
				curData = *cur_data_pos;
				pred_sz = cur_data_pos[-1] + cur_data_pos[-dim1_offset]+ cur_data_pos[-dim0_offset] - cur_data_pos[-dim1_offset - 1] - cur_data_pos[-dim0_offset - 1] - cur_data_pos[-dim0_offset - dim1_offset] + cur_data_pos[-dim0_offset - dim1_offset - 1];
				pred_reg = reg_params_pos[0] * i + reg_params_pos[params_offset_b] * i + reg_params_pos[params_offset_c] * i + reg_params_pos[params_offset_d];							
				err_sz += fabs(pred_sz - curData);
				err_reg += fabs(pred_reg - curData);

				bmi = block_size - i;
				cur_data_pos = data_pos + i*dim0_offset + i*dim1_offset + bmi;
				curData = *cur_data_pos;
				pred_sz = cur_data_pos[-1] + cur_data_pos[-dim1_offset]+ cur_data_pos[-dim0_offset] - cur_data_pos[-dim1_offset - 1] - cur_data_pos[-dim0_offset - 1] - cur_data_pos[-dim0_offset - dim1_offset] + cur_data_pos[-dim0_offset - dim1_offset - 1];
				pred_reg = reg_params_pos[0] * i + reg_params_pos[params_offset_b] * i + reg_params_pos[params_offset_c] * bmi + reg_params_pos[params_offset_d];							
				err_sz += fabs(pred_sz - curData);
				err_reg += fabs(pred_reg - curData);								

				cur_data_pos = data_pos + i*dim0_offset + bmi*dim1_offset + i;
				curData = *cur_data_pos;
				pred_sz = cur_data_pos[-1] + cur_data_pos[-dim1_offset]+ cur_data_pos[-dim0_offset] - cur_data_pos[-dim1_offset - 1] - cur_data_pos[-dim0_offset - 1] - cur_data_pos[-dim0_offset - dim1_offset] + cur_data_pos[-dim0_offset - dim1_offset - 1];
				pred_reg = reg_params_pos[0] * i + reg_params_pos[params_offset_b] * bmi + reg_params_pos[params_offset_c] * i + reg_params_pos[params_offset_d];							
				err_sz += fabs(pred_sz - curData);
				err_reg += fabs(pred_reg - curData);								

				cur_data_pos = data_pos + i*dim0_offset + bmi*dim1_offset + bmi;
				curData = *cur_data_pos;
				pred_sz = cur_data_pos[-1] + cur_data_pos[-dim1_offset]+ cur_data_pos[-dim0_offset] - cur_data_pos[-dim1_offset - 1] - cur_data_pos[-dim0_offset - 1] - cur_data_pos[-dim0_offset - dim1_offset] + cur_data_pos[-dim0_offset - dim1_offset - 1];
				pred_reg = reg_params_pos[0] * i + reg_params_pos[params_offset_b] * bmi + reg_params_pos[params_offset_c] * bmi + reg_params_pos[params_offset_d];							
				err_sz += fabs(pred_sz - curData);
				err_reg += fabs(pred_reg - curData);
			}
			reg_pred_err[pred_err_index] = err_reg;
			sz_pred_err[pred_err_index] = err_sz;
			pred_err_index ++;
		}
		reg_params_pos ++;
		sample_data_pos += block_size;
	}
}

void sz_sample_compress(float * sample_data, size_t num_sample_blocks, int block_size, double abs_eb, int use_mean, float mean, int quantization_intervals, float * reg_pred_err, float * sz_pred_err, float * reg_params, double * mse, float * bit_rate){

	size_t dim0_offset = num_sample_blocks * block_size * block_size;
	size_t dim1_offset = num_sample_blocks * block_size;
	double realPrecision = abs_eb;
	float noise = realPrecision * 1.12;
	size_t pred_err_index = 0;
	float * pred_buffer = (float *) malloc((block_size+1)*(block_size+1)*(block_size+1)*sizeof(float));
	float * pred_buffer_pos = NULL;
	float * block_data_pos_x = NULL;
	float * block_data_pos_y = NULL;
	float * block_data_pos_z = NULL;
	memset(pred_buffer, 0, (block_size+1)*(block_size+1)*(block_size+1)*sizeof(float));
	int pred_buffer_block_size = block_size + 1;
	int strip_dim0_offset = pred_buffer_block_size * pred_buffer_block_size;
	int strip_dim1_offset = pred_buffer_block_size;

    int * combined_sz_type = (int *) malloc(num_sample_blocks * block_size * block_size * block_size * sizeof(int));
    float * combined_sz_unpred_data = (float *) malloc (num_sample_blocks * block_size * block_size * block_size * sizeof(int));
	float rel_param_err = 0.025;
	double precision[4];
	precision[0] = rel_param_err * realPrecision / block_size;
	precision[1] = rel_param_err * realPrecision / block_size;
	precision[2] = rel_param_err * realPrecision / block_size;
	precision[3] = rel_param_err * realPrecision;
	float * reg_params_pos = reg_params;
	float last_coeffcients[4] = {0.0};
	int coeff_intvCapacity_sz = 65536;
	int coeff_intvRadius = 32768;
	int * coeff_type[4];
	int * coeff_result_type = (int *) malloc(num_sample_blocks*4*sizeof(int));
	float * coeff_unpred_data[4];
	float * coeff_unpredictable_data = (float *) malloc(num_sample_blocks*4*sizeof(float));
	for(int i=0; i<4; i++){
		coeff_type[i] = coeff_result_type + i * num_sample_blocks;
		coeff_unpred_data[i] = coeff_unpredictable_data + i * num_sample_blocks;
	}
	int coeff_index = 0;
	unsigned int coeff_unpredictable_count[4] = {0};

	SZ_Init(NULL);
	int * type = combined_sz_type;
    size_t combined_sz_unpred_count = 0;
	updateQuantizationInfo(quantization_intervals);
	int intvCapacity = exe_params->intvCapacity;
	int intvRadius = exe_params->intvRadius;	
	int intvCapacity_sz = intvCapacity - 2;

	size_t reg_count = 0;
	double final_mse = 0;
	float * sample_data_pos = sample_data;
	for(size_t i=0; i<num_sample_blocks; i++){
		// compress by sz
		// add 1 in x, y, z offset
		pred_buffer_pos = pred_buffer + pred_buffer_block_size*pred_buffer_block_size + pred_buffer_block_size + 1;
		block_data_pos_x = sample_data_pos;
		for(int ii=0; ii<block_size; ii++){
			block_data_pos_y = block_data_pos_x;
			for(int jj=0; jj<block_size; jj++){
				block_data_pos_z = block_data_pos_y;
				for(int kk=0; kk<block_size; kk++){
					*pred_buffer_pos = *block_data_pos_z;
					block_data_pos_z ++;
					pred_buffer_pos ++;
				}
				// add 1 in z offset
				pred_buffer_pos ++;
				block_data_pos_y += dim1_offset;
			}
			// add 1 in y offset
			pred_buffer_pos += pred_buffer_block_size;
			block_data_pos_x += dim0_offset;
		}

		if(reg_pred_err[pred_err_index] < sz_pred_err[pred_err_index] + noise * (block_size - 1) * 4){
			// regression
			{
				/*predict coefficients in current block via previous reg_block*/
				float cur_coeff;
				double diff, itvNum;
				for(int e=0; e<4; e++){
					cur_coeff = reg_params_pos[e*num_sample_blocks];
					diff = cur_coeff - last_coeffcients[e];
					itvNum = fabs(diff)/precision[e] + 1;
					if (itvNum < coeff_intvCapacity_sz){
						if (diff < 0) itvNum = -itvNum;
						coeff_type[e][coeff_index] = (int) (itvNum/2) + coeff_intvRadius;
						last_coeffcients[e] = last_coeffcients[e] + 2 * (coeff_type[e][coeff_index] - coeff_intvRadius) * precision[e];
						//ganrantee comporession error against the case of machine-epsilon
						if(fabs(cur_coeff - last_coeffcients[e])>precision[e]){	
							coeff_type[e][coeff_index] = 0;
							last_coeffcients[e] = cur_coeff;	
							coeff_unpred_data[e][coeff_unpredictable_count[e] ++] = cur_coeff;
						}					
					}
					else{
						coeff_type[e][coeff_index] = 0;
						last_coeffcients[e] = cur_coeff;
						coeff_unpred_data[e][coeff_unpredictable_count[e] ++] = cur_coeff;
					}
				}
				coeff_index ++;
			}
			float curData;
			float pred;
			double itvNum;
			double diff;
			size_t index = 0;
			float * cur_data_pos = pred_buffer + pred_buffer_block_size*pred_buffer_block_size + pred_buffer_block_size + 1;
			for(size_t ii=0; ii<block_size; ii++){
				for(size_t jj=0; jj<block_size; jj++){
					for(size_t kk=0; kk<block_size; kk++){
						curData = *cur_data_pos;
						pred = last_coeffcients[0] * ii + last_coeffcients[1] * jj + last_coeffcients[2] * kk + last_coeffcients[3];									
						diff = curData - pred;
						itvNum = fabs(diff)/realPrecision + 1;
						if (itvNum < intvCapacity){
							if (diff < 0) itvNum = -itvNum;
							type[index] = (int) (itvNum/2) + intvRadius;
							pred = pred + 2 * (type[index] - intvRadius) * realPrecision;
							//ganrantee comporession error against the case of machine-epsilon
							if(fabs(curData - pred)>realPrecision){	
								type[index] = 0;
								pred = curData;
								combined_sz_unpred_data[combined_sz_unpred_count ++] = curData;
							}		
						}
						else{
							type[index] = 0;
							pred = curData;
							combined_sz_unpred_data[combined_sz_unpred_count ++] = curData;
						}
						final_mse += (pred - curData) * (pred - curData);
						index ++;	
						cur_data_pos ++;
					}
					cur_data_pos ++;
				}
				cur_data_pos += pred_buffer_block_size;
			}
			reg_count ++;
		}
		else{
			// SZ
			if(use_mean){
				float * cur_data_pos = pred_buffer + pred_buffer_block_size*pred_buffer_block_size + pred_buffer_block_size + 1;
				float curData;
				float pred3D;
				double itvNum, diff;
				size_t index = 0;
				for(size_t ii=0; ii<block_size; ii++){
					for(size_t jj=0; jj<block_size; jj++){
						for(size_t kk=0; kk<block_size; kk++){

							curData = *cur_data_pos;
							if(fabs(curData - mean) <= realPrecision){
								type[index] = 1;
								*cur_data_pos = mean;
							}
							else
							{
								pred3D = cur_data_pos[-1] + cur_data_pos[-strip_dim1_offset]+ cur_data_pos[-strip_dim0_offset] - cur_data_pos[-strip_dim1_offset - 1]
										 - cur_data_pos[-strip_dim0_offset - 1] - cur_data_pos[-strip_dim0_offset - strip_dim1_offset] + cur_data_pos[-strip_dim0_offset - strip_dim1_offset - 1];
								diff = curData - pred3D;
								itvNum = fabs(diff)/realPrecision + 1;
								if (itvNum < intvCapacity_sz){
									if (diff < 0) itvNum = -itvNum;
									type[index] = (int) (itvNum/2) + intvRadius;
									*cur_data_pos = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
									//ganrantee comporession error against the case of machine-epsilon
									if(fabs(curData - *cur_data_pos)>realPrecision){	
										type[index] = 0;
										*cur_data_pos = curData;	
										combined_sz_unpred_data[combined_sz_unpred_count ++] = curData;
									}					
								}
								else{
									type[index] = 0;
									*cur_data_pos = curData;
									combined_sz_unpred_data[combined_sz_unpred_count ++] = curData;
								}
							}
							final_mse += (*cur_data_pos - curData) * (*cur_data_pos - curData);
							index ++;
							cur_data_pos ++;
						}
						cur_data_pos ++;
					}
					cur_data_pos += pred_buffer_block_size;
				}
			}
			else{
				float * cur_data_pos = pred_buffer + pred_buffer_block_size*pred_buffer_block_size + pred_buffer_block_size + 1;
				float curData;
				float pred3D;
				double itvNum, diff;
				size_t index = 0;
				for(size_t ii=0; ii<block_size; ii++){
					for(size_t jj=0; jj<block_size; jj++){
						for(size_t kk=0; kk<block_size; kk++){
							curData = *cur_data_pos;
							pred3D = cur_data_pos[-1] + cur_data_pos[-strip_dim1_offset]+ cur_data_pos[-strip_dim0_offset] - cur_data_pos[-strip_dim1_offset - 1]
									 - cur_data_pos[-strip_dim0_offset - 1] - cur_data_pos[-strip_dim0_offset - strip_dim1_offset] + cur_data_pos[-strip_dim0_offset - strip_dim1_offset - 1];
							diff = curData - pred3D;
							itvNum = fabs(diff)/realPrecision + 1;
							if (itvNum < intvCapacity_sz){
								if (diff < 0) itvNum = -itvNum;
								type[index] = (int) (itvNum/2) + intvRadius;
								*cur_data_pos = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
								//ganrantee comporession error against the case of machine-epsilon
								if(fabs(curData - *cur_data_pos)>realPrecision){	
									type[index] = 0;
									*cur_data_pos = curData;	
									combined_sz_unpred_data[combined_sz_unpred_count ++] = curData;
								}					
							}
							else{
								type[index] = 0;
								*cur_data_pos = curData;
								combined_sz_unpred_data[combined_sz_unpred_count ++] = curData;
							}											
							final_mse += (*cur_data_pos - curData) * (*cur_data_pos - curData);
							index ++;
							cur_data_pos ++;
						}
						cur_data_pos ++;
					}
					cur_data_pos += pred_buffer_block_size;
				}
			}// end use mean
		}// end else lorenzo
		type += block_size * block_size * block_size;
		pred_err_index ++;
		reg_params_pos ++;
		sample_data_pos += block_size;
	}
	free(pred_buffer);
	unsigned char * comp_data = (unsigned char *) malloc(num_sample_blocks*block_size*block_size*block_size*sizeof(float));
	unsigned char * comp_data_pos = comp_data;

	intToBytes_bigEndian(comp_data_pos, quantization_intervals);
	comp_data_pos += sizeof(int);
	memcpy(comp_data_pos, &use_mean, sizeof(unsigned char));
	comp_data_pos += sizeof(unsigned char);
	if(use_mean){
		memcpy(comp_data_pos, &mean, sizeof(float));
		comp_data_pos += sizeof(float);
	}
	if(reg_count > 0){
		for(int e=0; e<4; e++){
			int stateNum = 2*coeff_intvCapacity_sz;
			HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
			size_t nodeCount = 0;
			init(huffmanTree, coeff_type[e], reg_count);
			size_t i = 0;
			for (i = 0; i < huffmanTree->stateNum; i++)
				if (huffmanTree->code[i]) nodeCount++; 
			nodeCount = nodeCount*2-1;
			unsigned char *treeBytes;
			unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(huffmanTree, nodeCount, &treeBytes);
			intToBytes_bigEndian(comp_data_pos, treeByteSize);
			comp_data_pos += sizeof(int);
			intToBytes_bigEndian(comp_data_pos, nodeCount);
			comp_data_pos += sizeof(int);
			memcpy(comp_data_pos, treeBytes, treeByteSize);		
			comp_data_pos += treeByteSize;
			free(treeBytes);
			size_t typeArray_size = 0;
		
			encode(huffmanTree, coeff_type[e], reg_count, comp_data_pos + sizeof(size_t), &typeArray_size);
			sizeToBytes(comp_data_pos, typeArray_size);
			comp_data_pos += sizeof(size_t) + typeArray_size;
			intToBytes_bigEndian(comp_data_pos, coeff_unpredictable_count[e]);
			comp_data_pos += sizeof(int);
			memcpy(comp_data_pos, coeff_unpred_data[e], coeff_unpredictable_count[e]*sizeof(float));
			comp_data_pos += coeff_unpredictable_count[e]*sizeof(float);
			SZ_ReleaseHuffman(huffmanTree);
		}
	}
	*((size_t *) comp_data_pos) = combined_sz_unpred_count;
	comp_data_pos += sizeof(size_t);
	memcpy(comp_data_pos, combined_sz_unpred_data, combined_sz_unpred_count*sizeof(float));
	comp_data_pos += combined_sz_unpred_count*sizeof(float);
	{
		int stateNum = 2*quantization_intervals;
		HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
		size_t nodeCount = 0;
		init(huffmanTree, combined_sz_type, num_sample_blocks*block_size*block_size*block_size);
		size_t i = 0;
		for (i = 0; i < huffmanTree->stateNum; i++)
			if (huffmanTree->code[i]) nodeCount++; 
		nodeCount = nodeCount*2-1;
		unsigned char *treeBytes;
		unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(huffmanTree, nodeCount, &treeBytes);
		intToBytes_bigEndian(comp_data_pos, treeByteSize);
		comp_data_pos += sizeof(int);
		intToBytes_bigEndian(comp_data_pos, nodeCount);
		comp_data_pos += sizeof(int);
		memcpy(comp_data_pos, treeBytes, treeByteSize);
		comp_data_pos += treeByteSize;
		free(treeBytes);
		size_t typeArray_size = 0;
		encode(huffmanTree, combined_sz_type, num_sample_blocks*block_size*block_size*block_size, comp_data_pos + sizeof(size_t), &typeArray_size);
		*((size_t *) comp_data_pos) = typeArray_size;
		comp_data_pos += sizeof(size_t);
		comp_data_pos += typeArray_size;
		SZ_ReleaseHuffman(huffmanTree);
	}
	unsigned char * compressed_comp_data = NULL;
	size_t totalEncodeSize = comp_data_pos - comp_data;
	size_t lossless_outsize = sz_lossless_compress(confparams_cpr->losslessCompressor, confparams_cpr->gzipMode, comp_data, totalEncodeSize, &compressed_comp_data);
	SZ_Finalize();
	*bit_rate = lossless_outsize * 32.0 / (num_sample_blocks * block_size * block_size * block_size * sizeof(float));
	*mse = final_mse/(num_sample_blocks * block_size * block_size * block_size);
	free(comp_data);
	free(compressed_comp_data);
	if(coeff_result_type) free(coeff_result_type);
	if(coeff_unpredictable_data) free(coeff_unpredictable_data);
	if(combined_sz_type) free(combined_sz_type);
	if(combined_sz_unpred_data) free(combined_sz_unpred_data);
}

float * get_sample_data(float * data, size_t r1, size_t r2, size_t r3, int block_size, int sample_distance, size_t *sample_blocks){
	size_t num_x = r1 / block_size;
	size_t num_y = r2 / block_size;
	size_t num_z = r3 / block_size;
	int sample_block_size = MIN(num_y, num_z);
	sample_block_size = MIN(sample_block_size, num_x);
	int sample_nx = num_x / sample_block_size;
	int sample_ny = num_y / sample_block_size;
	int sample_nz = num_z / sample_block_size;

	sample_nx = (sample_nx - 1) / sample_distance + 1;
	sample_ny = (sample_ny - 1) / sample_distance + 1;
	sample_nz = (sample_nz - 1) / sample_distance + 1;
	size_t num_sample_blocks = sample_nx * sample_ny * sample_nz * sample_block_size * 4;
	*sample_blocks = num_sample_blocks;

	float * sample_data = (float *) malloc(num_sample_blocks*block_size*block_size*block_size*sizeof(float));
	float * sample_data_pos = sample_data;
	size_t sample_data_dim0_offset = num_sample_blocks * block_size * block_size;
	size_t sample_data_dim1_offset = num_sample_blocks * block_size;
	size_t * offset_x = (size_t *) malloc(sample_block_size*sizeof(size_t));
	size_t * offset_y = (size_t *) malloc(sample_block_size*sizeof(size_t));
	size_t * offset_z = (size_t *) malloc(sample_block_size*sizeof(size_t));
	for(int i=0; i<sample_block_size; i++){
		offset_x[i] = i*block_size*r2*r3;
		offset_y[i] = i*block_size*r3;
		offset_z[i] = i*block_size;
	}
	// printf("%d %d %d, %d %zu\n", sample_nx, sample_ny, sample_nz, sample_block_size, num_sample_blocks);
	for(size_t ix=0; ix<sample_nx; ix+=sample_distance){
		for(size_t iy=0; iy<sample_ny; iy+=sample_distance){
			for(size_t iz=0; iz<sample_nz; iz+=sample_distance){
				float * block_data = data + ix*block_size*sample_block_size*r2*r3 + iy*block_size*sample_block_size*r3 + iz*block_size*sample_block_size;
				// sample inside blocks, sample the four diagonals
				// sample the first diagonal [i, i, i]
				for(int i=0; i<sample_block_size; i++){
					float * cur_block_data_pos = block_data + offset_x[i] + offset_y[i] + offset_z[i];
					float * cur_sample_data_pos = sample_data_pos;
					for(int ii=0; ii<block_size; ii++){
						for(int jj=0; jj<block_size; jj++){
							for(int kk=0; kk<block_size; kk++){
								*(cur_sample_data_pos ++) = *(cur_block_data_pos ++);
							}
							cur_sample_data_pos += sample_data_dim1_offset - block_size;
							cur_block_data_pos += r3 - block_size;
						}
						cur_sample_data_pos += sample_data_dim0_offset - sample_data_dim1_offset*block_size;
						cur_block_data_pos += r2*r3 - r3*block_size;
					}
					sample_data_pos += block_size;
				}
				// sample the second diagonal [i, i, n-i]
				for(int i=0; i<sample_block_size; i++){
					float * cur_block_data_pos = block_data + offset_x[i] + offset_y[i] + offset_z[sample_block_size - 1 - i];
					float * cur_sample_data_pos = sample_data_pos;
					for(int ii=0; ii<block_size; ii++){
						for(int jj=0; jj<block_size; jj++){
							for(int kk=0; kk<block_size; kk++){
								*(cur_sample_data_pos ++) = *(cur_block_data_pos ++);
							}
							cur_sample_data_pos += sample_data_dim1_offset - block_size;
							cur_block_data_pos += r3 - block_size;
						}
						cur_sample_data_pos += sample_data_dim0_offset - sample_data_dim1_offset*block_size;
						cur_block_data_pos += r2*r3 - r3*block_size;
					}
					sample_data_pos += block_size;
				}
				// sample the thrid diagonal [i, n-i, i]
				for(int i=0; i<sample_block_size; i++){
					float * cur_block_data_pos = block_data + offset_x[i] + offset_y[sample_block_size - 1 - i] + offset_z[i];
					float * cur_sample_data_pos = sample_data_pos;
					for(int ii=0; ii<block_size; ii++){
						for(int jj=0; jj<block_size; jj++){
							for(int kk=0; kk<block_size; kk++){
								*(cur_sample_data_pos ++) = *(cur_block_data_pos ++);
							}
							cur_sample_data_pos += sample_data_dim1_offset - block_size;
							cur_block_data_pos += r3 - block_size;
						}
						cur_sample_data_pos += sample_data_dim0_offset - sample_data_dim1_offset*block_size;
						cur_block_data_pos += r2*r3 - r3*block_size;
					}
					sample_data_pos += block_size;
				}
				// sample the last diagonal [i, n-i, n-i]
				for(int i=0; i<sample_block_size; i++){
					float * cur_block_data_pos = block_data + offset_x[i] + offset_y[sample_block_size - 1 - i] + offset_z[sample_block_size - 1 - i];
					float * cur_sample_data_pos = sample_data_pos;
					for(int ii=0; ii<block_size; ii++){
						for(int jj=0; jj<block_size; jj++){
							for(int kk=0; kk<block_size; kk++){
								*(cur_sample_data_pos ++) = *(cur_block_data_pos ++);
							}
							cur_sample_data_pos += sample_data_dim1_offset - block_size;
							cur_block_data_pos += r3 - block_size;
						}
						cur_sample_data_pos += sample_data_dim0_offset - sample_data_dim1_offset*block_size;
						cur_block_data_pos += r2*r3 - r3*block_size;
					}
					sample_data_pos += block_size;
				}
			}
		}
	}
	// printf("sample done\n");
	free(offset_x);
	free(offset_y);
	free(offset_z);
	return sample_data;
}

signed int * transform_sample_data(float * sample_data, size_t num_sample_blocks, int block_size, int * block_exp){
	signed int * ordered_fp = (signed int *) malloc(num_sample_blocks * block_size * block_size * block_size * sizeof(signed int));
	size_t dim0_offset = num_sample_blocks * block_size * block_size;
	size_t dim1_offset = num_sample_blocks * block_size;
	float * sample_data_pos = sample_data;
	const int bs = 4;
	signed int * ordered_fp_pos = ordered_fp;
	int * block_exp_pos = block_exp;
	float block_data[64];
	signed int block_data_fp[64];
	int wlen = 32;
	int intbits = 2;
	init_blockid();
	for(size_t i=0; i<num_sample_blocks; i++){
		// for each sample, transform 2x2x2 4x4x4 blocks
		for(int i=0; i<2; i++){
			for(int j=0; j<2; j++){
				for(int k=0; k<2; k++){
					float * data_pos = sample_data_pos + i*bs*dim0_offset + j*bs*dim1_offset + k*bs;
					float max_ele = 0;
					float * block_data_pos = block_data;
					for(int ii=0; ii<4; ii++){
						for(int jj=0; jj<4; jj++){
							for(int kk=0; kk<4; kk++){
								float cur_data = *(data_pos + ii*dim0_offset + jj*dim1_offset + kk);
								*(block_data_pos++) = cur_data;
								max_ele = MAX(max_ele, fabs(cur_data));
							}
						}
					}
					// to fixed point
					int e;
				    frexp(max_ele, &e);
				    *(block_exp_pos ++) = e;
					{
						for(int ii=0; ii<64; ii++){
					    	float cur_data = ldexp(block_data[ii], wlen - intbits - e);
					    	block_data_fp[ii] = (long int) cur_data;
					    }
					}
					// transform
					// transform along z
					for(int ii=0; ii<4; ii++){
						for(int jj=0; jj<4; jj++){
							fp_fwd_lift(block_data_fp + ii*16 + 4*jj, 1);
						}
					}
					// transform along y
					for(int ii=0; ii<4; ii++){
						for(int kk=0; kk<4; kk++){
							fp_fwd_lift(block_data_fp + ii*16 + kk, 4);
						}
					}
					// transform along x
					for(int jj=0; jj<4; jj++){
						for(int kk=0; kk<4; kk++){
							fp_fwd_lift(block_data_fp + jj*4 + kk, 16);
						}
					}
					// reorder
					signed int * block_data_fp_pos = block_data_fp;
					for(int ii=0; ii<4; ii++){
						for(int jj=0; jj<4; jj++){
							for(int kk=0; kk<4; kk++){
								int id = block_id[ii][jj][kk];
								ordered_fp_pos[id] = *(block_data_fp_pos ++);
							}
						}
					}
					ordered_fp_pos += 64;
				}
			}
		}
		sample_data_pos += block_size;
	}
	return ordered_fp;
}

float * zfp_decompress(unsigned char * comp_data, size_t comp_data_size, size_t comp_data_size_before_lossless, size_t r1, size_t r2, size_t r3){
	size_t num_x, num_y, num_z;
	int block_size = 8;
	num_x = (r1 - 1) / block_size + 1;
	num_y = (r2 - 1) / block_size + 1;
	num_z = (r3 - 1) / block_size + 1;
	size_t num_blocks = num_x * num_y * num_z;
	size_t num_elements = r1 * r2 * r3;
	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;

	int * combined_zfp_exp = NULL;
	MemoryBitStream stream;
	int sz_compress_num = 0;
	float * combined_zfp_sz_compress_coeff_unpred_data = NULL;
	int * combined_zfp_type = NULL;
	int * combined_zfp_sz_compress_coeff_type = NULL;
	float * combined_zfp_unpred_data = NULL;
	int zfp_coeff_intv_capacity = 65536;
	int zfp_coeff_intv_radius = 32768;
	int zfp_wrap_intv_capacity = 65536;
	int zfp_wrap_intv_radius = 32768;
	int zfp_factor = 10;

	SZ_Init(NULL);
	unsigned char * comp_data_uncompressed = NULL;
	size_t lossless_output = sz_lossless_decompress(confparams_cpr->losslessCompressor, comp_data, comp_data_size, &comp_data_uncompressed, comp_data_size_before_lossless);

	unsigned char * comp_data_pos = comp_data_uncompressed;
	double realPrecision = *((double *) comp_data_pos);
	comp_data_pos += sizeof(double);

	zfp_factor = *((int *) comp_data_pos);
	comp_data_pos += sizeof(float);
	// block exp
	size_t compressed_exp_size = *((size_t *) comp_data_pos);
	comp_data_pos += sizeof(size_t);
	combined_zfp_exp = (int *)SZ_decompress(SZ_INT32, comp_data_pos, compressed_exp_size, 0, 0, 0, 0, num_blocks * 8);
	comp_data_pos += compressed_exp_size;
	// embedded encoding
	size_t embedded_coding_size = *((size_t *) comp_data_pos);
	comp_data_pos += sizeof(size_t);
	stream.open(comp_data_pos, embedded_coding_size);
	comp_data_pos += embedded_coding_size;		
	// sz_compress num
	sz_compress_num = *((int *) comp_data_pos);
	comp_data_pos += sizeof(int);
	if(sz_compress_num){
		// unpred data
		size_t combined_zfp_sz_compress_coeff_unpred_count = *((size_t *) comp_data_pos);
		comp_data_pos += sizeof(size_t);
		combined_zfp_sz_compress_coeff_unpred_data = (float *) comp_data_pos;
		comp_data_pos += combined_zfp_sz_compress_coeff_unpred_count * sizeof(float);

		int stateNum = zfp_coeff_intv_capacity;
		HuffmanTree* huffmanTree = createHuffmanTree(stateNum);	
		int tree_size = bytesToInt_bigEndian(comp_data_pos);
		comp_data_pos += sizeof(int);
		int nodeCount = bytesToInt_bigEndian(comp_data_pos);
		comp_data_pos += sizeof(int);
		node root = reconstruct_HuffTree_from_bytes_anyStates(huffmanTree, comp_data_pos, nodeCount);
		comp_data_pos += tree_size;
		size_t type_array_size = *((size_t *)comp_data_pos);
		comp_data_pos += sizeof(size_t);
		combined_zfp_sz_compress_coeff_type = (int *) malloc(sz_compress_num*num_blocks*8 * sizeof(int));
		decode(comp_data_pos, sz_compress_num*num_blocks*8, root, combined_zfp_sz_compress_coeff_type);
		comp_data_pos += type_array_size;
		SZ_ReleaseHuffman(huffmanTree);
	}
	// wrap
	{
		// type array
		int stateNum = zfp_wrap_intv_capacity;
		HuffmanTree* huffmanTree = createHuffmanTree(stateNum);	
		int tree_size = bytesToInt_bigEndian(comp_data_pos);
		comp_data_pos += sizeof(int);
		int nodeCount = bytesToInt_bigEndian(comp_data_pos);
		comp_data_pos += sizeof(int);
		node root = reconstruct_HuffTree_from_bytes_anyStates(huffmanTree, comp_data_pos, nodeCount);
		comp_data_pos += tree_size;

		size_t type_array_size = *((size_t *)comp_data_pos);
		comp_data_pos += sizeof(size_t);
		combined_zfp_type = (int *) malloc(num_blocks*block_size*block_size*block_size * sizeof(int));
		decode(comp_data_pos, num_blocks*block_size*block_size*block_size, root, combined_zfp_type);
		comp_data_pos += type_array_size;
		SZ_ReleaseHuffman(huffmanTree);
		// unpred
		size_t zfp_wrap_unpred_count = *((size_t *) comp_data_pos);
		comp_data_pos += sizeof(size_t);
		combined_zfp_unpred_data = (float *) comp_data_pos;
		comp_data_pos += zfp_wrap_unpred_count * sizeof(float);
	}
	// decompress
    ptrdiff_t zfp_offsets[8];
    for(int i=0; i<2; i++){
    	for(int j=0; j<2; j++){
    		for(int k=0; k<2; k++){
    			zfp_offsets[4*i+2*j+k] = i*4*dim0_offset + j*4*dim1_offset + k*4;
    		}
    	}
    }
    // int zfp_order[8] = {0, 1, 2, 4, 3, 5, 6, 7};
    int zfp_order[8] = {6, 2, 0, 4, 5, 1, 3, 7};
    float * zfp_data_pos;
    float * zfp_last_sz_compressed_coeff = (float *) malloc(sz_compress_num*sizeof(float));
    memset(zfp_last_sz_compressed_coeff, 0, sz_compress_num * sizeof(float));
    int * block_exp_pos = combined_zfp_exp;

	int emin = INT_MIN;
	double zfp_precision = realPrecision * zfp_factor;
	if (zfp_precision > 0) {
		frexp(zfp_precision, &emin);
		emin--;
	}
	double zfp_coeff_eb = ldexp(1.0, emin) / 64; // 4^3 = 64 is for the maximum loss in transform

	int gp_element_num[9] = {1, 3, 6, 10, 12, 12, 10, 6, 4};
	int gp_element_offset[9] = {0, 1, 4, 10, 20, 32, 44, 54, 60};

	int64 count = 0x46acca631ull;
	int sz_compress_group = 0;
	{
		int tmp = sz_compress_num;
		while(tmp > 0){
			tmp -= gp_element_num[sz_compress_group];
			count >>= 4;
			sz_compress_group ++;
		}
	}
	int wlen = 32;
	int intbits = 2;

	signed int block_data_fp[64];
	signed int ordered_fp[64];
	int * sz_compress_coeff_type = combined_zfp_sz_compress_coeff_type;
	int * zfp_wrap_type = combined_zfp_type;
	size_t zfp_sz_compress_coeff_index = 0;
	size_t zfp_wrap_unpred_count = 0;
	size_t combined_zfp_sz_compress_coeff_unpred_count = 0;
	float * dec_data = (float *) malloc(num_elements * sizeof(float));
	float * data_pos;
	init_blockid();
    for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				data_pos = dec_data + i*block_size * dim0_offset + j*block_size * dim1_offset + k*block_size;
				// printf("%d %d %d\n", i, j, k);
				// mzfp
				for(int z=0; z<8; z++){
					// coefficient analysis
					{
						int e = *block_exp_pos;
						// restore sz_compressed coeff
						for(int num=0; num<sz_compress_num; num++){
							// decompress sz compressed coeff
							float cur_dec_data;
							if(sz_compress_coeff_type[zfp_sz_compress_coeff_index] == 0){
								cur_dec_data = combined_zfp_sz_compress_coeff_unpred_data[combined_zfp_sz_compress_coeff_unpred_count ++];
							}
							else{
								cur_dec_data = zfp_last_sz_compressed_coeff[num] + 2 * (sz_compress_coeff_type[zfp_sz_compress_coeff_index] - zfp_coeff_intv_radius) * zfp_coeff_eb;
							}
							zfp_sz_compress_coeff_index ++;
							zfp_last_sz_compressed_coeff[num] = cur_dec_data;
					    	float cur_data = ldexp(cur_dec_data, wlen - intbits - e);
					    	ordered_fp[num] = (long int) cur_data;
						}
						// embedded decoding
						int maxprec = MAX(0, e - emin + 8);
						embedded_decoding(stream, ordered_fp + sz_compress_num, maxprec, count, 64 - sz_compress_num);
						// reorder
						for(int ii=0; ii<4; ii++){
							for(int jj=0; jj<4; jj++){
								for(int kk=0; kk<4; kk++){
									int id = block_id[ii][jj][kk];
									*(block_data_fp + ii*16 + jj*4 + kk) = ordered_fp[id];
								}
							}
						}
					}				
					// transform along x
					for(int jj=0; jj<4; jj++){
						for(int kk=0; kk<4; kk++){
							fp_inv_lift(block_data_fp + jj*4 + kk, 16);
						}
					}
					// transform along y
					for(int ii=0; ii<4; ii++){
						for(int kk=0; kk<4; kk++){
							fp_inv_lift(block_data_fp + ii*16 + kk, 4);
						}
					}
					// transform along z
					for(int ii=0; ii<4; ii++){
						for(int jj=0; jj<4; jj++){
							fp_inv_lift(block_data_fp + ii*16 + 4*jj, 1);
						}
					}
					// restore
					float block_data[64];
					{
						int e = *block_exp_pos;					
						for(int ii=0; ii<64; ii++){
					    	signed int fp = block_data_fp[ii];
					    	block_data[ii] = ldexp((float)fp, intbits - wlen + e);
					    }
					}

					zfp_data_pos = data_pos + zfp_offsets[zfp_order[z]];
					{
						int wrap_index = 0;
						for(int ii=0; ii<4; ii++){
							for(int jj=0; jj<4; jj++){
								for(int kk=0; kk<4; kk++){
									float cur_data = *(block_data + ii*16 + jj*4 + kk);
									// quantization
									if(zfp_wrap_type[wrap_index] == 0){
										cur_data = combined_zfp_unpred_data[zfp_wrap_unpred_count ++];
									}
									else{
										cur_data = cur_data + 2 * (zfp_wrap_type[wrap_index] - zfp_wrap_intv_radius) * realPrecision;
									}
									wrap_index ++;
									*(zfp_data_pos + ii*r2*r3 + jj*r3 + kk) = cur_data;
								}
							}
						}
						zfp_wrap_type += 64;
					}
					
				    block_exp_pos ++;
				}
			}
		}
	}
	if(zfp_last_sz_compressed_coeff) free(zfp_last_sz_compressed_coeff);
	if(combined_zfp_exp) free(combined_zfp_exp);
	if(combined_zfp_type) free(combined_zfp_type);
	if(combined_zfp_sz_compress_coeff_type) free(combined_zfp_sz_compress_coeff_type);
	if(comp_data_uncompressed) free(comp_data_uncompressed);
	SZ_Finalize();
	return dec_data;
}

float * zfp_decompress_blocksize_4(unsigned char * comp_data, size_t comp_data_size, size_t comp_data_size_before_lossless, size_t r1, size_t r2, size_t r3){
	size_t num_x, num_y, num_z;
	int block_size = 4;
	num_x = (r1 - 1) / block_size + 1;
	num_y = (r2 - 1) / block_size + 1;
	num_z = (r3 - 1) / block_size + 1;
	size_t num_blocks = num_x * num_y * num_z;
	size_t num_elements = r1 * r2 * r3;
	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;

	int * combined_zfp_exp = NULL;
	MemoryBitStream stream;
	int sz_compress_num = 0;
	float * combined_zfp_sz_compress_coeff_unpred_data = NULL;
	int * combined_zfp_type = NULL;
	int * combined_zfp_sz_compress_coeff_type = NULL;
	float * combined_zfp_unpred_data = NULL;
	int zfp_coeff_intv_capacity = 65536;
	int zfp_coeff_intv_radius = 32768;
	int zfp_wrap_intv_capacity = 65536;
	int zfp_wrap_intv_radius = 32768;
	int zfp_factor = 10;

	SZ_Init(NULL);
	unsigned char * comp_data_uncompressed = NULL;
	size_t lossless_output = sz_lossless_decompress(confparams_cpr->losslessCompressor, comp_data, comp_data_size, &comp_data_uncompressed, comp_data_size_before_lossless);

	unsigned char * comp_data_pos = comp_data_uncompressed;
	double realPrecision = *((double *) comp_data_pos);
	comp_data_pos += sizeof(double);

	zfp_factor = *((int *) comp_data_pos);
	comp_data_pos += sizeof(float);
	// block exp
	size_t compressed_exp_size = *((size_t *) comp_data_pos);
	comp_data_pos += sizeof(size_t);
	combined_zfp_exp = (int *)SZ_decompress(SZ_INT32, comp_data_pos, compressed_exp_size, 0, 0, 0, 0, num_blocks);
	comp_data_pos += compressed_exp_size;
	// embedded encoding
	size_t embedded_coding_size = *((size_t *) comp_data_pos);
	comp_data_pos += sizeof(size_t);
	stream.open(comp_data_pos, embedded_coding_size);
	comp_data_pos += embedded_coding_size;		
	// sz_compress num
	sz_compress_num = *((int *) comp_data_pos);
	comp_data_pos += sizeof(int);
	if(sz_compress_num){
		// unpred data
		size_t combined_zfp_sz_compress_coeff_unpred_count = *((size_t *) comp_data_pos);
		comp_data_pos += sizeof(size_t);
		combined_zfp_sz_compress_coeff_unpred_data = (float *) comp_data_pos;
		comp_data_pos += combined_zfp_sz_compress_coeff_unpred_count * sizeof(float);

		int stateNum = zfp_coeff_intv_capacity;
		HuffmanTree* huffmanTree = createHuffmanTree(stateNum);	
		int tree_size = bytesToInt_bigEndian(comp_data_pos);
		comp_data_pos += sizeof(int);
		int nodeCount = bytesToInt_bigEndian(comp_data_pos);
		comp_data_pos += sizeof(int);
		node root = reconstruct_HuffTree_from_bytes_anyStates(huffmanTree, comp_data_pos, nodeCount);
		comp_data_pos += tree_size;
		size_t type_array_size = *((size_t *)comp_data_pos);
		comp_data_pos += sizeof(size_t);
		combined_zfp_sz_compress_coeff_type = (int *) malloc(sz_compress_num*num_blocks* sizeof(int));
		decode(comp_data_pos, sz_compress_num*num_blocks, root, combined_zfp_sz_compress_coeff_type);
		comp_data_pos += type_array_size;
		SZ_ReleaseHuffman(huffmanTree);
	}
	// wrap
	{
		// type array
		int stateNum = zfp_wrap_intv_capacity;
		HuffmanTree* huffmanTree = createHuffmanTree(stateNum);	
		int tree_size = bytesToInt_bigEndian(comp_data_pos);
		comp_data_pos += sizeof(int);
		int nodeCount = bytesToInt_bigEndian(comp_data_pos);
		comp_data_pos += sizeof(int);
		node root = reconstruct_HuffTree_from_bytes_anyStates(huffmanTree, comp_data_pos, nodeCount);
		comp_data_pos += tree_size;

		size_t type_array_size = *((size_t *)comp_data_pos);
		comp_data_pos += sizeof(size_t);
		combined_zfp_type = (int *) malloc(num_blocks*block_size*block_size*block_size * sizeof(int));
		decode(comp_data_pos, num_blocks*block_size*block_size*block_size, root, combined_zfp_type);
		comp_data_pos += type_array_size;
		SZ_ReleaseHuffman(huffmanTree);
		// unpred
		size_t zfp_wrap_unpred_count = *((size_t *) comp_data_pos);
		comp_data_pos += sizeof(size_t);
		combined_zfp_unpred_data = (float *) comp_data_pos;
		comp_data_pos += zfp_wrap_unpred_count * sizeof(float);
	}
	// decompress
    float * zfp_last_sz_compressed_coeff = (float *) malloc(sz_compress_num*sizeof(float));
    memset(zfp_last_sz_compressed_coeff, 0, sz_compress_num * sizeof(float));
    int * block_exp_pos = combined_zfp_exp;

	int emin = INT_MIN;
	double zfp_precision = realPrecision * zfp_factor;
	if (zfp_precision > 0) {
		frexp(zfp_precision, &emin);
		emin--;
	}
	double zfp_coeff_eb = ldexp(1.0, emin) / 64; // 4^3 = 64 is for the maximum loss in transform

	int gp_element_num[9] = {1, 3, 6, 10, 12, 12, 10, 6, 4};
	int gp_element_offset[9] = {0, 1, 4, 10, 20, 32, 44, 54, 60};

	int64 count = 0x46acca631ull;
	int sz_compress_group = 0;
	{
		int tmp = sz_compress_num;
		while(tmp > 0){
			tmp -= gp_element_num[sz_compress_group];
			count >>= 4;
			sz_compress_group ++;
		}
	}
	int wlen = 32;
	int intbits = 2;

	signed int block_data_fp[64];
	signed int ordered_fp[64];
	int * sz_compress_coeff_type = combined_zfp_sz_compress_coeff_type;
	int * zfp_wrap_type = combined_zfp_type;
	size_t zfp_sz_compress_coeff_index = 0;
	size_t zfp_wrap_unpred_count = 0;
	size_t combined_zfp_sz_compress_coeff_unpred_count = 0;
	float * dec_data = (float *) malloc(num_elements * sizeof(float));
	float * data_pos;
	init_blockid();
    for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				data_pos = dec_data + i*block_size * dim0_offset + j*block_size * dim1_offset + k*block_size;
				// printf("%d %d %d\n", i, j, k);
				// mzfp
				// coefficient analysis
				{
					int e = *block_exp_pos;
					// restore sz_compressed coeff
					for(int num=0; num<sz_compress_num; num++){
						// decompress sz compressed coeff
						float cur_dec_data;
						if(sz_compress_coeff_type[zfp_sz_compress_coeff_index] == 0){
							cur_dec_data = combined_zfp_sz_compress_coeff_unpred_data[combined_zfp_sz_compress_coeff_unpred_count ++];
						}
						else{
							cur_dec_data = zfp_last_sz_compressed_coeff[num] + 2 * (sz_compress_coeff_type[zfp_sz_compress_coeff_index] - zfp_coeff_intv_radius) * zfp_coeff_eb;
						}
						zfp_sz_compress_coeff_index ++;
						zfp_last_sz_compressed_coeff[num] = cur_dec_data;
				    	float cur_data = ldexp(cur_dec_data, wlen - intbits - e);
				    	ordered_fp[num] = (long int) cur_data;
					}
					// embedded decoding
					int maxprec = MAX(0, e - emin + 8);
					embedded_decoding(stream, ordered_fp + sz_compress_num, maxprec, count, 64 - sz_compress_num);
					// reorder
					for(int ii=0; ii<4; ii++){
						for(int jj=0; jj<4; jj++){
							for(int kk=0; kk<4; kk++){
								int id = block_id[ii][jj][kk];
								*(block_data_fp + ii*16 + jj*4 + kk) = ordered_fp[id];
							}
						}
					}
				}				
				// transform along x
				for(int jj=0; jj<4; jj++){
					for(int kk=0; kk<4; kk++){
						fp_inv_lift(block_data_fp + jj*4 + kk, 16);
					}
				}
				// transform along y
				for(int ii=0; ii<4; ii++){
					for(int kk=0; kk<4; kk++){
						fp_inv_lift(block_data_fp + ii*16 + kk, 4);
					}
				}
				// transform along z
				for(int ii=0; ii<4; ii++){
					for(int jj=0; jj<4; jj++){
						fp_inv_lift(block_data_fp + ii*16 + 4*jj, 1);
					}
				}
				// restore
				float block_data[64];
				{
					int e = *block_exp_pos;					
					for(int ii=0; ii<64; ii++){
				    	signed int fp = block_data_fp[ii];
				    	block_data[ii] = ldexp((float)fp, intbits - wlen + e);
				    }
				}

				float * zfp_data_pos = data_pos;
				{
					int wrap_index = 0;
					for(int ii=0; ii<4; ii++){
						for(int jj=0; jj<4; jj++){
							for(int kk=0; kk<4; kk++){
								float cur_data = *(block_data + ii*16 + jj*4 + kk);
								// quantization
								if(zfp_wrap_type[wrap_index] == 0){
									cur_data = combined_zfp_unpred_data[zfp_wrap_unpred_count ++];
								}
								else{
									cur_data = cur_data + 2 * (zfp_wrap_type[wrap_index] - zfp_wrap_intv_radius) * realPrecision;
								}
								wrap_index ++;
								*(zfp_data_pos + ii*r2*r3 + jj*r3 + kk) = cur_data;
							}
						}
					}
					zfp_wrap_type += 64;
				}
				
			    block_exp_pos ++;
			}
		}
	}
	if(zfp_last_sz_compressed_coeff) free(zfp_last_sz_compressed_coeff);
	if(combined_zfp_exp) free(combined_zfp_exp);
	if(combined_zfp_type) free(combined_zfp_type);
	if(combined_zfp_sz_compress_coeff_type) free(combined_zfp_sz_compress_coeff_type);
	if(comp_data_uncompressed) free(comp_data_uncompressed);
	SZ_Finalize();
	return dec_data;
}

unsigned char * compress_block(float * data, size_t r3, size_t r2, size_t r1, double eb, size_t * out_size, size_t * out_size_before_lossless, int * select){
	size_t num_elements = r3*r2*r1;
	float value_range = 0;
    {
	    float max_ele = data[0], min_ele = data[0];
	    for(int i=0; i<num_elements; i++){
			max_ele = MAX(max_ele, data[i]);
			min_ele = MIN(min_ele, data[i]);
	    }
	    value_range = max_ele - min_ele;
	}
	double realPrecision = value_range * eb;
	// int sz_compress_num = estimate_zfp_sz_compress_num(data, r1, r2, r3, value_range * eb, 16);
	// estimate zfp nrmse & cr
	int zfp_factor = 1;
	double nrmse = 0;
	float bit_rate = 0;
	double mse[ZFP_FACTOR_SAMPLE_NUM];
	int sz_compress_num = 0;
	int sample_sz_compress_num[ZFP_FACTOR_SAMPLE_NUM];
	float bit_rates[ZFP_FACTOR_SAMPLE_NUM];
	// get sample data
	size_t num_sample_blocks = 0;
	int block_size = 8;
	int sample_distance = 1;
	float * sample_data = get_sample_data(data, r1, r2, r3, block_size, sample_distance, &num_sample_blocks);
	// block exp for zfp sample, *8 for 8 blocks
	int * block_exp = (int *) malloc(num_sample_blocks * 8 * sizeof(int));
	signed int * ordered_fp = transform_sample_data(sample_data, num_sample_blocks, block_size, block_exp);
	int best_pos = 0;
	float best_ratio = 0;
	for(int i=0; i<ZFP_FACTOR_SAMPLE_NUM; i++){
		zfp_sample_compress(sample_data, ordered_fp, block_exp, num_sample_blocks, value_range * eb, zfp_factor, &sample_sz_compress_num[i], &mse[i], &bit_rates[i]);
		if(i){
			if((bit_rates[i] >= bit_rates[i-1]) && (mse[i] >= mse[i-1])) break;
			float bit_rate_decay = log2(bit_rates[0] / bit_rates[i]);
			float mse_increase = log2(mse[i] / mse[0]);
			float cur_ratio = bit_rate_decay / mse_increase;
			if(cur_ratio > best_ratio){
				best_pos = i;
				best_ratio = cur_ratio;
			}
		}
		zfp_factor *= 2;
	}
	int best_zfp_factor = 1 << best_pos;
	float zfp_psnr = 20*log10(value_range)-10*log10(mse[best_pos]);
	float zfp_bit_rate = bit_rates[best_pos];
	sz_compress_num = sample_sz_compress_num[best_pos];
	// printf("Current best pos %d, eb %d\n", best_pos, best_zfp_factor);
	// printf("ZFP Bit_rate = %.4g\tNRMSE = %.4g\tPSNR = %.4g\n", zfp_bit_rate, sqrt(best_mse)/value_range, zfp_psnr);

	// SZ statistics
	SZ_Init(NULL);
	unsigned int quantization_intervals;
	unsigned char use_mean = 0;
	float dense_pos;
	if(exe_params->optQuantMode==1)
	{
		float sz_sample_correct_freq = -1;//0.5; //-1
		float mean_flush_freq;
		quantization_intervals = optimize_intervals_float_3D_with_freq_and_dense_pos(data, r1, r2, r3, realPrecision, &dense_pos, &sz_sample_correct_freq, &mean_flush_freq);
		if(mean_flush_freq > 0.5 || mean_flush_freq > sz_sample_correct_freq) use_mean = 1;
		updateQuantizationInfo(quantization_intervals);
	}	
	else{
		quantization_intervals = exe_params->intvCapacity;
	}
	float mean = 0;
	if(use_mean){
		// compute mean
		double sum = 0.0;
		size_t mean_count = 0;
		for(size_t i=0; i<num_elements; i++){
			if(fabs(data[i] - dense_pos) < realPrecision){
				sum += data[i];
				mean_count ++;
			}
		}
		if(mean_count > 0) mean = sum / mean_count;
	}
	SZ_Finalize();
	float * reg_pred_err = (float *) malloc(num_sample_blocks * sizeof(float));
	float * sz_pred_err = (float *) malloc(num_sample_blocks * sizeof(float));
	float * reg_params = (float *) malloc(num_sample_blocks * 4 * sizeof(float));
	sz_compute_sample_reg_coefficient(sample_data, num_sample_blocks, block_size, use_mean, mean, reg_params, reg_pred_err, sz_pred_err);
	double cur_eb = eb;
	double last_eb = eb;
	float cur_psnr;
	float cur_bit_rate;
	double cur_mse;
	sz_sample_compress(sample_data, num_sample_blocks, block_size, value_range * cur_eb, use_mean, mean, quantization_intervals, reg_pred_err, sz_pred_err, reg_params, &cur_mse, &cur_bit_rate);
	cur_psnr = 20*log10(value_range)-10*log10(cur_mse);
	float last_psnr = cur_psnr;
	float last_bit_rate = cur_bit_rate;
	bool use_sz = false;
	bool use_zfp = false;
	if(cur_psnr < zfp_psnr){
		while(cur_psnr < zfp_psnr){
			if(cur_bit_rate > zfp_bit_rate){
				use_zfp = true;
				break;
			}
			last_psnr = cur_psnr;
			last_bit_rate = cur_bit_rate;
			last_eb = cur_eb;
			cur_eb /= 2;
			sz_sample_compress(sample_data, num_sample_blocks, block_size, value_range * cur_eb, use_mean, mean, quantization_intervals, reg_pred_err, sz_pred_err, reg_params, &cur_mse, &cur_bit_rate);
			cur_psnr = 20*log10(value_range)-10*log10(cur_mse);
		}
	}
	else{
		while(cur_psnr > zfp_psnr){
			if(cur_bit_rate < zfp_bit_rate){
				use_sz = true;
				break;
			}
			last_psnr = cur_psnr;
			last_bit_rate = cur_bit_rate;
			last_eb = cur_eb;
			cur_eb *= 2;
			sz_sample_compress(data, num_sample_blocks, block_size, value_range * cur_eb, use_mean, mean, quantization_intervals, reg_pred_err, sz_pred_err, reg_params, &cur_mse, &cur_bit_rate);
			cur_psnr = 20*log10(value_range)-10*log10(cur_mse);
		}
	}
	// printf("SZ eb %.4g: Bit_rate = %.4g\tPSNR = %.4g\n", last_eb, last_bit_rate, last_psnr);
	// printf("SZ eb %.4g: Bit_rate = %.4g\tPSNR = %.4g\n", cur_eb, cur_bit_rate, cur_psnr);
	free(reg_pred_err);
	free(sz_pred_err);
	free(reg_params);
	free(sample_data);
	free(block_exp);
	free(ordered_fp);

	float interpolated_sz_psnr = (cur_psnr - last_psnr) / (cur_bit_rate - last_bit_rate) * (zfp_bit_rate - last_bit_rate) + last_psnr; 
	unsigned char * comp_data = NULL;
	//use_sz = true;use_zfp = false;
	//use_zfp = true;
	if(use_sz) interpolated_sz_psnr = zfp_psnr + 1;
	if(use_zfp) interpolated_sz_psnr = zfp_psnr - 1;
	if(interpolated_sz_psnr > zfp_psnr){
		SZ_Init(NULL);
		comp_data = SZ_compress_args(SZ_FLOAT, (void *) data, out_size, REL, 0, eb, 0, 0, 0, r1, r2, r3);
		SZ_Finalize();
		*select = 0;
	}
	else{
		// printf("sz_compress_num: %d\n", sz_compress_num);
		comp_data = zfp_compress_block_size_4(data, r1, r2, r3, value_range*eb, best_zfp_factor, sz_compress_num, out_size, out_size_before_lossless, &mse[0], &bit_rates[0]);
		// comp_data = zfp_compress(data, r1, r2, r3, value_range*eb, best_zfp_factor, sz_compress_num, out_size, out_size_before_lossless, &mse[0], &bit_rates[0]);
		*select = 1;
	}
	return comp_data;
}

float * decompress_block(unsigned char * comp_data, size_t comp_data_size, size_t comp_data_size_before_lossless, int select, int r3, int r2, int r1){
	if(select) return zfp_decompress_blocksize_4(comp_data, comp_data_size, comp_data_size_before_lossless, r1, r2, r3);
	// if(select) return zfp_decompress(comp_data, comp_data_size, comp_data_size_before_lossless, r1, r2, r3);
	else{
		SZ_Init(NULL);
		float * dec_data = (float *) SZ_decompress(SZ_FLOAT, comp_data, comp_data_size, 0, 0, r1, r2, r3);
		SZ_Finalize();
		return dec_data;
	}
}




