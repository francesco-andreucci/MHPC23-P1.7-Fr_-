# Molecular Dynamics Simulation Code

This package is designed for simulating atoms using a simplified MD (Molecular Dynamics) code with multi-threading parallelization, specifically optimized for a Lennard-Jones potential.

## Project Overview

The `examples` directory includes three sets of example input decks, and the `reference` directory contains the corresponding outputs.

## Compilation Instructions

### Prerequisites
In the CMakeList.txt make sure that the following options are configured based on your requirements:

1. **option(USE_OMP "Building with OPENMP" ON/OFF)**: Activate OpenMP support.
2. **option(USE_MPI "Building with MPI" ON/OFF)**: Activate MPI support.
3. **option(ENABLE_TESTING "Enable building unit tests" ON)/OFF**: Enable to compile tests.

### Compilation Steps
This package contains simplified MD code with multi-threading parallelization for simulating atoms with a Lennard-Jones potential. The examples directory contains 3 sets of example input decks and the reference directory the corresponding outputs. 
How to compile:

1. Navigate to the project's root folder.
2. Run the following commands:

```bash
   cmake -S . -B build
   cmake --build build
```
or use the bmake.sh script (./bmake.sh build). To run the code, cd to the created build folder then you can run:

1. Serial one: 
```bash
    ./main < input_file.inp
```
2. OpenMP support
```bash
    OMP_NUM_THREADS=#threads ./main < input_file.inp
```
3. MPI support
```bash
    mpirun -np #mpitaks ./main < input_file.inp
```
4. Hybrid OpenMP & MPI
```bash
    OMP_NUM_THREADS=#threads mpirun -np #mpitaks ./main < input_file.inp
```
Where #threads, and  #mpitaks stands for the desired number of threads and the number of MPI tasks, respectively.

We would like to mention that this project uses gtest library. Unit tests are configured using googletest, and they automatically detect whether OpenMP and/or MPI are present or not and test accordingly. All the test can be run by typing "make test" from the build folder.

NB. On macOS with Clang, additional settings are applied to address specific warnings and OpenMP support.

Our report is in the file Benchmark_report_MHPC_FZB-1.pdf.

We split the work in the following way:
1. Serial optimization: Behzad Salmassian
2. MPI implementation: Francesco Andreucci
3. OpenMP implementation: Zakaria Dahbi
4. Hybrid MPI-OpenMP: Francesco and Zakaria

The various commits begin with the name of the commiter. The commits up to Sunday 29/11 (included) refer to the group assignments, the ones after refer to the individual tasks.