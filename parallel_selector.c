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
#include "sz.h"
#include "rw.h"
#include "sz_zfp_selector.h"
#include "mpi.h"

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
	
	SZ_Init(NULL);

	if (world_rank == 0) printf ("Start parallel compressing ... \n");
	if (world_rank == 0) printf("size: %d\n", world_size);
	double start, end;
	double costReadOri = 0.0, costReadZip = 0.0, costWriteZip = 0.0, costWriteOut = 0.0, costComp = 0.0, costDecomp = 0.0;

	MPI_Barrier(MPI_COMM_WORLD);

	int num_vars = atoi(argv[2]);
	// NYX
	int nyx_num_vars = 6;
	char nyx_file[6][50] ={"velocity_x.dat", "velocity_y.dat", "velocity_z.dat", "dark_matter_density_log10.dat", "temperature.dat", "baryon_density_log10.dat"};
	double nyx_rel_bound[6] = {5e-3, 5e-3, 5e-3, 9e-2, 1e-2, 5e-2};
	// Hurricane
	int hurricane_num_vars = 13;
	char hurricane_file[13][50] ={"Uf48_truncated.bin.dat", "Vf48_truncated.bin.dat", "Wf48_truncated.bin.dat", "TCf48_truncated.bin.dat", "Pf48_truncated.bin.dat", "QVAPORf48_truncated.bin.dat", "CLOUDf48_log10_truncated.bin.dat", "QCLOUDf48_log10_truncated.bin.dat", "QICEf48_log10_truncated.bin.dat", "QRAINf48_log10_truncated.bin.dat", "QSNOWf48_log10_truncated.bin.dat", "QGRAUPf48_log10_truncated.bin.dat", "PRECIPf48_log10_truncated.bin.dat"};
	double hurricane_rel_bound[13] = {6e-3, 3e-3, 3e-3, 4e-3, 4e-3, 5e-3, 3.8e-2, 4e-2, 4e-2, 2.5e-2, 2.7e-2, 2.6e-2, 3e-2};
	// SCALE 
	int scale_num_vars = 12;
	char scale_file[12][50] ={"U_truncated.bin.dat", "V_truncated.bin.dat", "W_truncated.bin.dat", "T_truncated.bin.dat", "RH_truncated.bin.dat", "PRES_truncated.bin.dat", "QC_log10_truncated.bin.dat", "QR_log10_truncated.bin.dat", "QI_log10_truncated.bin.dat", "QS_log10_truncated.bin.dat", "QG_log10_truncated.bin.dat", "QV_log10_truncated.bin.dat"};
	double scale_rel_bound[12] = {4e-3, 8e-3, 2e-2, 4e-3, 5e-3, 5e-3, 2.5e-2, 5.6e-2, 2e-2, 8e-2, 3e-2, 3e-1};

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
		SZ_Finalize();
		MPI_Finalize();
		return 0;
	}	
	size_t compressed_size[13];
	int select[13];
	size_t comp_data_size_before_lossless[13];

	char folder[50] = "/lcrc/project/ECP-EZ/public/compression/datasets";
	char filename[100];
	char zip_filename[100];
	// char out_filename[100];
	size_t inSize, outSize; 
	size_t nbEle;
	int status;
	float * dataIn;

	size_t est_compressed_size = r1 * r2 * r3 * sizeof(float) * num_vars / 30;
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
		if (world_rank == 0) printf ("Compressing %s\n", filename);
		MPI_Barrier(MPI_COMM_WORLD);
		if(world_rank == 0) start = MPI_Wtime();
		unsigned char *bytesOut = compress_block(dataIn, r1, r2, r3, rel_bound[i], &compressed_size[i], &comp_data_size_before_lossless[i], &select[i]); 
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

	sprintf(zip_filename, "%s/%d/selector.out", folder, folder_index);	// Write Compressed Data
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
		float *dataOut = decompress_block(compressed_output_pos, compressed_size[i], comp_data_size_before_lossless[i], select[i], r1, r2, r3);
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
		printf ("Selector Finish parallel compressing, total compression ratio %.4g.\n", 1.0*r1*r2*r3*sizeof(float)*num_vars / total_size);
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

	SZ_Finalize();

	MPI_Finalize();

	return 0;
}
