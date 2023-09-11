# PDC Summer School
Repository with codes, exercises, and the project for the PDC Summer School 2023.

The repo will be continuously updated during the school, so make sure to pull any new changes. 

For questions about a particular exercise session, write in slack or ask us in person during the hands-on.

Have a great summer school!

https://coderefinery.github.io/research-software-engineering/

To clone this repo:
`git clone --recurse-submodules https://github.com/KTH-HPC/pdc-summer-school.git`

To pull any changes (get new updates, changes and exercises that are added):
`git pull https://github.com/KTH-HPC/pdc-summer-school.git`

To also load the content of the linked repositories (such as the GPU exercises), use the following git command:
`git submodule update --init --recursive`

# Running the project

To run the OpenMP version, use
export OMP_NUM_THREADS=<num_threads>
`./energy_storms_omp <size> <storm_1_file> [ <storm_i_file>] ...`

To run the MPI version, use
`mpirun -n <num_threads> ./energy_storms_mpi <size> <storm_1_file> [ <storm_i_file>] ...`

To run the CUDA version, use
`./energy_storms_cuda <num_threads> <size> <storm_1_file> [ <storm_i_file>] ...`

# Benchmarks

To benchmark the versions on the provided test cases, use
`./run_tests.sh <num_threads> <file_with_test_cases>`
(Python 3 required)

e.g
`./run_tests.sh 8 all_tests.txt`
