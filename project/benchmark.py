#!/usr/bin/env python

import argparse
import os

from collections import namedtuple

RunResults = namedtuple('RunResults', ['time', 'results'])

versions = {"seq": "Sequential", "omp": "OpenMP", "mpi": "MPI", "cuda": "CUDA"}

def createparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', "-r", help="Directory with benchmark results")
    return parser

def parse_max_string(max_string):
    res_spl = max_string.strip().split()
    result_list = [(int(res_spl[2 * i + 1]), float(res_spl[2 * i + 2])) for i in range((len(res_spl) - 1) // 2)]
    return result_list

def cmp_to_baseline(baseline, parallel):
    err_threshold = 0.1
    # print(baseline, parallel)
    if len(baseline) != len(parallel):
        return False
    # assert(len(baseline) % 2 == 0)
    for i in range(len(baseline)):
        if baseline[i][0] != parallel[i][0]:
            return False
        if abs(baseline[i][1] - parallel[i][1]) > err_threshold:
            return False
    return True

def validate_results(test_to_results):
    failed_tests = {}
    for test, results in test_to_results.items():
        baseline = results["seq"]
        for version, run_results in results.items():
            if version != "seq":
                if version not in failed_tests:
                    failed_tests[version] = 0
                if not cmp_to_baseline(baseline.results, run_results.results):
                    failed_tests[version] += 1
    for version, num in failed_tests.items():
        print("{} failed tests for {}".format(num, versions[version]))

def time_benchmark(test_to_results):
    print("Test number\t{}\n".format("\t".join([versions[version] for version in versions.keys()])))
    for test in sorted(test_to_results.keys()):
        results = test_to_results[test]
        time_arr = [str(test)]
        times = [str(results[version].time) for version in versions.keys()]
        time_str = "\t".join(time_arr + times)
        print(time_str)

def benchmark_results(test_to_results):
    validate_results(test_to_results)
    time_benchmark(test_to_results)

def parse_results(results_dir):
    test_to_run_results = {}
    for path in os.listdir(results_dir):
        if path.endswith(".txt"):
            version, num = path.split(".")[0].split("_")
            time = None
            max_list = None
            with open(os.path.join(results_dir, path), "r") as in_handle:
                for line in in_handle:
                    if line.startswith("Time"):
                        time = float(line.strip().split()[1])
                    if line.startswith("Result"):
                        max_list = parse_max_string(line)
            run_results = RunResults(time=time, results=max_list)
            if num not in test_to_run_results:
                test_to_run_results[num] = {}
            test_to_run_results[num][version] = run_results
    return test_to_run_results


parser = createparser()
args = parser.parse_args()

test_to_results = parse_results(args.results)
benchmark_results(test_to_results)
