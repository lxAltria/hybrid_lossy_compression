#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstddef>
#include "sz.h"
#include "rw.h"
#include "sz_zfp_selector.h"
#include <iostream>

#define MAX(a, b) a>b?a:b
#define MIN(a, b) a<b?a:b

int main(int argc, char ** argv){
    int status;
    size_t nbEle;
    float * data = readFloatData(argv[1], &nbEle, &status);
    int n3 = atoi(argv[2]);
    int n2 = atoi(argv[3]);
    int n1 = atoi(argv[4]);
    double eb = atof(argv[5]);
    printf("start compress\n");
    size_t out_size;
    size_t comp_data_size_before_lossless;
    int select;
    unsigned char * comp_data = compress_block(data, n3, n2, n1, eb, &out_size, &comp_data_size_before_lossless, &select);
    std::string compressors[2] = {"SZ", "ZFP"};
    printf("\nSelect=%d %s, compressed_size: %zu, ratio: %.4g\n\n", select, compressors[select].c_str(), out_size, nbEle*sizeof(float)*1.0 / out_size);
    writeByteData(comp_data, out_size, strcat(argv[1], ".select"), &status);
    float * dec_data = decompress_block(comp_data, out_size, comp_data_size_before_lossless, select, n3, n2, n1);
    writeFloatData_inBytes(dec_data, n1*n2*n3, strcat(argv[1], ".out"), &status);
	printf("Decompress all over\n");
    // verify
	{
        float * ori_data = data;
        float * data = dec_data;
        size_t i = 0;
        float Max = 0, Min = 0, diffMax = 0;
        Max = ori_data[0];
        Min = ori_data[0];
        diffMax = fabs(data[0] - ori_data[0]);
        size_t k = 0;
        double sum1 = 0, sum2 = 0;
        for (i = 0; i < nbEle; i++)
        {
            sum1 += ori_data[i];
            sum2 += data[i];
        }
        double mean1 = sum1/nbEle;
        double mean2 = sum2/nbEle;

        double sum3 = 0, sum4 = 0;
        double sum = 0, prodSum = 0, relerr = 0;

        double maxpw_relerr = 0; 
        for (i = 0; i < nbEle; i++)
        {
            if (Max < ori_data[i]) Max = ori_data[i];
            if (Min > ori_data[i]) Min = ori_data[i];
            
            float err = fabs(data[i] - ori_data[i]);
            if(ori_data[i]!=0 && fabs(ori_data[i])>1)
            {
                relerr = err/fabs(ori_data[i]);
                if(maxpw_relerr<relerr)
                    maxpw_relerr = relerr;
            }

            if (diffMax < err)
                diffMax = err;
            prodSum += (ori_data[i]-mean1)*(data[i]-mean2);
            sum3 += (ori_data[i] - mean1)*(ori_data[i]-mean1);
            sum4 += (data[i] - mean2)*(data[i]-mean2);
            sum += err*err; 
        }
        double std1 = sqrt(sum3/nbEle);
        double std2 = sqrt(sum4/nbEle);
        double ee = prodSum/nbEle;
        double acEff = ee/std1/std2;

        double mse = sum/nbEle;
        double range = Max - Min;
        double psnr = 20*log10(range)-10*log10(mse);
        double nrmse = sqrt(mse)/range;
        // double compressionRatio = 1.0*(n1 * n2 * n3) * sizeof(float) / compressed_size;
        double compressionRatio = 1.0*(n1 * n2 * n3) * sizeof(float) / out_size;//1.0*32 / (32 - total_aver_tz - total_aver_lz);

        printf ("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
        printf ("Max absolute error = %.10f\n", diffMax);
        printf ("Max relative error = %f\n", diffMax/(Max-Min));
        printf ("Max pw relative error = %f\n", maxpw_relerr);
        printf ("PSNR = %f, NRMSE= %.20G\n", psnr,nrmse);
        printf ("acEff=%f\n", acEff);   
        printf ("compressionRatio=%f\n", compressionRatio);
    }
    free(comp_data);
    free(data);
    free(dec_data);
}
