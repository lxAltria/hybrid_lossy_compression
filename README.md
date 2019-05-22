Source codes for a hybrid lossy compression framework

Instrunctions on our solution:
1. Download SZ (https://github.com/disheng222/SZ.git). It is better to use the version committed on Oct. 2018 to ensure correctness.
2. Find "confparams_cpr->predThreshold = 0.99" in SZ/sz/src/conf.c and change it to "confparams_cpr->predThreshold = 0.999".
3. Find "block_size = 6" and in SZ/sz/src/sz_float.c and change all the occurence to "block_size = 8". 
4. Install SZ
5. Install ZFP-0.3.1
6. Compile the code (link to the libraries built above) using the shell script.
Optional
7. We add a block-wise version to improve memory overhead in the rebuttal. To compile this version, do "cp sz_zfp_selector_block_version.cpp sz_zfp_selector.cpp" and re-compile the code. Please note that this only works for the sequential code for now ("sz_zfp_select.cpp", we also only evaluate the compression ratio instead of parallel performance in the paper).

Instructions on comparison:
1. Download and install SZ-2.0 (http://dstats.net/download/http://www.mcs.anl.gov/~shdi/download/sz-2.0.0.0.tar.gz)
2. Download and install ZFP-0.5.4 (https://computation.llnl.gov/projects/floating-point-compression/download/zfp-0.5.4.tar.gz)
3. Compile the corresponding code ("parallel_sz_2.0.c" and "")

Preprocessing on dataset:
1. Download the datasets from https://sdrbench.github.io. Please download the data with logarithmic transform (except SCALE as it only has original data) as this kind of data is for visualization purpose.
2. Run the processing script for Hurricane and SCALE. As SCALE data with logarithmic transform is not available, we do it by ourself in the script. 
3. Run the experiments: sequential: 	./select velocity_x.dat 512 512 512 1e-3 			#NYX
					./select CLOUDf48_log10_truncated.bin.dat 496 496 96 1e-3 	#Hurricane
					./select T_truncated.bin.dat 1200 1200 96 1e-3 			#SCALE
			parallel:	mpirun -np 32 ./parallel_selector sz.config 6 512 512 512	#NYX
					mpirun -np 32 ./parallel_selector sz.config 13 496 496 96	#Hurricane
					mpirun -np 32 ./parallel_selector sz.config 12 1200 1200 96	#SCALE

