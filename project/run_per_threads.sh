#!/usr/bin/env bash

for threads in 2 4 8 12; 
do
	export OMP_NUM_THREADS=$threads
	echo $threads
	echo "OpenMP"
	./energy_storms_omp 1000000 test_files/test_07_a1M_p5k_w1 test_files/test_07_a1M_p5k_w2 test_files/test_07_a1M_p5k_w3 test_files/test_07_a1M_p5k_w4 | grep "Time"
	echo "MPI"
	mpirun -n $threads "./energy_storms_mpi" 1000000 test_files/test_07_a1M_p5k_w1 test_files/test_07_a1M_p5k_w2 test_files/test_07_a1M_p5k_w3 test_files/test_07_a1M_p5k_w4 < /dev/null | grep "Time" 
	echo "CUDA"
	./energy_storms_cuda $threads 1000000 test_files/test_07_a1M_p5k_w1 test_files/test_07_a1M_p5k_w2 test_files/test_07_a1M_p5k_w3 test_files/test_07_a1M_p5k_w4 | grep "Time"
done 
