#!/usr/bin/env bash

while read test_case; 
do
	./energy_storms_seq $test_case | sed '/^Time/d' > test_results/energy_storms_seq_result
	$1 $test_case | sed '/^Time/d' > test_results/parallel_result
	cmp test_results/parallel_result test_results/energy_storms_seq_result
done < "test_cmds.txt"
