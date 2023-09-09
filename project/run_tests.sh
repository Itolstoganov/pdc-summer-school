#!/usr/bin/env bash

test_num=1
export OMP_NUM_THREADS
while read test_case; 
do
	echo $test_case
	echo $test_num
	./energy_storms_seq $test_case > test_results/seq_$test_num.txt
	./energy_storms_omp $test_case > test_results/omp_$test_num.txt
	mpirun -n 9 "./energy_storms_mpi" $test_case > test_results/mpi_$test_num.txt < /dev/null
	let "test_num += 1"
done < $1

python3 benchmark.py --results test_results
