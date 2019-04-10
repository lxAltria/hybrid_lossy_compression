/**
 *  @file test_compress_parallel.c
 *  @author Dingwen Tao
 *  @date January, 2017
 *  @brief This is an example of using compression interface in parallel
 *  (C) 2017 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zfp.h"
#include "mpi.h"
#include "rw.h"

unsigned char * zfp_compress_3D(float * array, double tolerance, size_t r1, size_t r2, size_t r3, size_t *out_size){
	int status = 0;    /* return value: 0 = success */
	zfp_type type;     /* array scalar type */
	zfp_field* field;  /* array meta data */
	zfp_stream* zfp;   /* compressed stream */
	void* buffer;      /* storage for compressed stream */
	size_t bufsize;    /* byte size of compressed buffer */
	bitstream* stream; /* bit stream to write to or read from */
	size_t zfpsize;    /* byte size of compressed stream */

	/* allocate meta data for the 3D array a[nz][ny][nx] */
	type = zfp_type_float;
	field = zfp_field_3d(array, type, r3, r2, r1);

	/* allocate meta data for a compressed stream */
	zfp = zfp_stream_open(NULL);

	/* set compression mode and parameters via one of three functions */
	/*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
	/*  zfp_stream_set_precision(zfp, precision); */
	zfp_stream_set_accuracy(zfp, tolerance);

	/* allocate buffer for compressed data */
	bufsize = zfp_stream_maximum_size(zfp, field);
	buffer = malloc(bufsize);

	/* associate bit stream with allocated buffer */
	stream = stream_open(buffer, bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	zfp_stream_rewind(zfp);

	zfpsize = zfp_compress(zfp, field);
    if (!zfpsize) {
      fprintf(stderr, "compression failed\n");
      status = 1;
    }	

	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);
	*out_size = zfpsize;
	return (unsigned char *)buffer;
}

float * zfp_decompress_3D(unsigned char * comp_data, double tolerance, size_t buffer_size, size_t r1, size_t r2, size_t r3){
	int status = 0;    /* return value: 0 = success */
	zfp_type type;     /* array scalar type */
	zfp_field* field;  /* array meta data */
	zfp_stream* zfp;   /* compressed stream */
	void* buffer;      /* storage for compressed stream */
	size_t bufsize;    /* byte size of compressed buffer */
	bitstream* stream; /* bit stream to write to or read from */
	size_t zfpsize;    /* byte size of compressed stream */

	/* allocate meta data for the 3D array a[nz][ny][nx] */
	float * array = (float *) malloc(r1 * r2 * r3 * sizeof(float));
	type = zfp_type_float;
	field = zfp_field_3d(array, type, r3, r2, r1);

	/* allocate meta data for a compressed stream */
	zfp = zfp_stream_open(NULL);

	/* set compression mode and parameters via one of three functions */
	/*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
	/*  zfp_stream_set_precision(zfp, precision); */
	zfp_stream_set_accuracy(zfp, tolerance);

	/* allocate buffer for compressed data */
	bufsize = zfp_stream_maximum_size(zfp, field);
	// buffer = malloc(bufsize);
	buffer = (void *) comp_data;
	bufsize = buffer_size;

	/* associate bit stream with allocated buffer */
	stream = stream_open(buffer, bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	zfp_stream_rewind(zfp);

    if (!zfp_decompress(zfp, field)) {
      fprintf(stderr, "decompression failed\n");
      status = 1;
    }
	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);
	return array;
}

// USAGE
// mpirun -np 16 parallel sz.config folder_num r3 r2 r1
int main(int argc, char * argv[])
{

	size_t r5=0,r4=0,r3=0,r2=0,r1=0;
	char *cfgFile;

	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);	
	
	if(argc < 3)
	{
		printf("Test case: testfloat_compress [config_file] [srcFilePath] [dimension sizes...]\n");
		printf("Example: testfloat_compress sz.config testfloat_8_8_128.dat 8 8 128\n");
		exit(0);
	}

	cfgFile=argv[1];
	
	if(argc>=4)
	  r1 = atoi(argv[3]); //8
	if(argc>=5)
	  r2 = atoi(argv[4]); //8
	if(argc>=6)
	  r3 = atoi(argv[5]); //128
	if(argc>=7)
	  r4 = atoi(argv[6]);
	if(argc>=8)
	  r5 = atoi(argv[7]);
	
	if (world_rank == 0) printf ("Start parallel compressing ... \n");
	if (world_rank == 0) printf("size: %d\n", world_size);
	double start, end;
	double costReadOri = 0.0, costReadZip = 0.0, costWriteZip = 0.0, costWriteOut = 0.0, costComp = 0.0, costDecomp = 0.0;

	MPI_Barrier(MPI_COMM_WORLD);

	int num_vars = atoi(argv[2]);
	// NYX
	int nyx_num_vars = 6;
	char nyx_file[6][50] ={"velocity_x.dat", "velocity_y.dat", "velocity_z.dat", "dark_matter_density_log10.dat", "temperature.dat", "baryon_density_log10.dat"};
	double nyx_rel_bound[6] = {4e6, 4e6, 4e6, 20, 4e5, 8};
	// Hurricane
	int hurricane_num_vars = 13;
	char hurricane_file[13][50] ={"Uf48_truncated.bin.dat", "Vf48_truncated.bin.dat", "Wf48_truncated.bin.dat", "TCf48_truncated.bin.dat", "Pf48_truncated.bin.dat", "QVAPORf48_truncated.bin.dat", "CLOUDf48_log10_truncated.bin.dat", "QCLOUDf48_log10_truncated.bin.dat", "QICEf48_log10_truncated.bin.dat", "QRAINf48_log10_truncated.bin.dat", "QSNOWf48_log10_truncated.bin.dat", "QGRAUPf48_log10_truncated.bin.dat", "PRECIPf48_log10_truncated.bin.dat"};
	double hurricane_rel_bound[13] = {8, 4, 0.5, 4, 400, 1e-3, 4, 16, 8, 16, 16, 16, 16};
	// SCALE 
	int scale_num_vars = 12;
	char scale_file[12][50] ={"U_truncated.bin.dat", "V_truncated.bin.dat", "W_truncated.bin.dat", "T_truncated.bin.dat", "RH_truncated.bin.dat", "PRES_truncated.bin.dat", "QC_log10_truncated.bin.dat", "QR_log10_truncated.bin.dat", "QI_log10_truncated.bin.dat", "QS_log10_truncated.bin.dat", "QG_log10_truncated.bin.dat", "QV_log10_truncated.bin.dat"};
	double scale_rel_bound[12] = {4, 8, 4, 8, 16, 12000, 16, 8, 16, 16, 16, 8};

	// assignment
	char file[13][50];
	double * rel_bound;
	if(num_vars == nyx_num_vars){
		for(int i=0; i<num_vars; i++) strcpy(file[i], nyx_file[i]);
		rel_bound = nyx_rel_bound;
	}
	else if(num_vars == hurricane_num_vars){
		for(int i=0; i<num_vars; i++) strcpy(file[i], hurricane_file[i]);
		rel_bound = hurricane_rel_bound;
	}
	else if(num_vars == scale_num_vars){
		for(int i=0; i<num_vars; i++) strcpy(file[i], scale_file[i]);
		rel_bound = scale_rel_bound;
	}
	else{
		printf("No such variablem, exit\n");
		MPI_Finalize();
		return 0;
	}
	size_t compressed_size[13];


	char folder[50] = "/lcrc/project/ECP-EZ/public/compression/datasets";
	char filename[100];
	char zip_filename[100];
	// char out_filename[100];
	size_t inSize, outSize; 
	size_t nbEle;
	int status;
	float * dataIn;

	size_t est_compressed_size = r1 * r2 * r3 * sizeof(float) * num_vars / 25;
	unsigned char * compressed_output = (unsigned char *) malloc(est_compressed_size);
	unsigned char * compressed_output_pos = compressed_output;
	int folder_index = world_rank;
	for(int i=0; i<num_vars; i++){
		sprintf(filename, "%s/%d/%s", folder, folder_index, file[i]);
		sprintf(zip_filename, "%s/%d/%s.sz", folder, folder_index, file[i]);
		// Read Input Data
		if(world_rank == 0){
			start = MPI_Wtime();
			dataIn = readFloatData(filename, &nbEle, &status);
			end = MPI_Wtime();
			printf("file %s read time: %.2f\n", filename, end - start);
			start = MPI_Wtime();
			MPI_Bcast(&nbEle, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
			MPI_Bcast(dataIn, nbEle, MPI_FLOAT, 0, MPI_COMM_WORLD);
			end = MPI_Wtime();
			printf("broadcast time: %.2f\n", end - start);
		}
		else{
			MPI_Bcast(&nbEle, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
			dataIn = (float *) malloc(nbEle * sizeof(float));
			MPI_Bcast(dataIn, nbEle, MPI_FLOAT, 0, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(world_rank == 0){
			end = MPI_Wtime();
			costReadOri += end - start;
		}
		
		// Compress Input Data
		size_t out_size;
		if (world_rank == 0) printf ("Compressing %s\n", filename);
		MPI_Barrier(MPI_COMM_WORLD);
		if(world_rank == 0) start = MPI_Wtime();
		unsigned char * bytesOut = zfp_compress_3D(dataIn, rel_bound[i], r3, r2, r1, &compressed_size[i]);
		MPI_Barrier(MPI_COMM_WORLD);
		if(world_rank == 0){
			end = MPI_Wtime();
			costComp += end - start;
		}
		free (dataIn);
		memcpy(compressed_output_pos, bytesOut, compressed_size[i]);
		compressed_output_pos += compressed_size[i];
		free(bytesOut);
	}

	sprintf(zip_filename, "%s/%d/zfp.out", folder, folder_index);	// Write Compressed Data
	size_t total_size = compressed_output_pos - compressed_output;
	// Write Compressed Data
	MPI_Barrier(MPI_COMM_WORLD);
	if(world_rank == 0) start = MPI_Wtime();
	writeByteData(compressed_output, total_size, zip_filename, &status);
	MPI_Barrier(MPI_COMM_WORLD);
	if(world_rank == 0){
		end = MPI_Wtime();
		costWriteZip += end - start;
	}
	free(compressed_output);
	// Read Compressed Data
	MPI_Barrier(MPI_COMM_WORLD);
	if(world_rank == 0) start = MPI_Wtime();
	compressed_output = readByteData(zip_filename, &inSize, &status);
	MPI_Barrier(MPI_COMM_WORLD);
	if(world_rank == 0){
		end = MPI_Wtime();
		costReadZip += end - start;
	}
	compressed_output_pos = compressed_output;

	for(int i=0; i<num_vars; i++){
		// Decompress Compressed Data
		MPI_Barrier(MPI_COMM_WORLD);
		if(world_rank == 0) start = MPI_Wtime();
		float * dataOut = zfp_decompress_3D(compressed_output_pos, rel_bound[i], compressed_size[i], r3, r2, r1);
		MPI_Barrier(MPI_COMM_WORLD);
		if(world_rank == 0){
			end = MPI_Wtime();
			costDecomp += end - start; 
		}
		compressed_output_pos += compressed_size[i];
		free(dataOut);
	}
	free(compressed_output);

	if (world_rank == 0)
	{
		printf ("ZFP Finish parallel compressing, total compression ratio %.4g.\n", 1.0*r1*r2*r3*sizeof(float)*num_vars / total_size);
		printf("Separate ratios: ");
		for(int i=0; i<num_vars; i++){
			printf("%.4g ", 1.0*r1*r2*r3*sizeof(float) / compressed_size[i]);
		}
		printf("\n");
		printf ("Timecost of reading original files = %.2f seconds\n", costReadOri);
		printf ("Timecost of reading compressed files = %.2f seconds\n", costReadZip);
		printf ("Timecost of writing compressed files = %.2f seconds\n", costWriteZip);
		printf ("Timecost of writing decompressed files = %.2f seconds\n", costWriteOut);
		printf ("Timecost of compressing using %d processes = %.2f seconds\n", world_size, costComp);
		printf ("Timecost of decompressing using %d processes = %.2f seconds\n\n", world_size, costDecomp);
	}


	MPI_Finalize();

	return 0;
}
