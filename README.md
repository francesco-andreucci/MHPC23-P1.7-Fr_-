This package contains simplified MD code with multi-threading
parallelization for simulating atoms with a Lennard-Jones potential.

The examples directory contains 3 sets of example input decks
and the reference directory the corresponding outputs.

How to compile:

1) Turn on/off the relevant options:
   USE_OPENMP: activate OpenMP support
   USE_MPI: activate MPI support
   ENABLE_TESTING: to compile the tests
2) From the root folder: cmake -S . -B build
3) cmake --build build

Unit tests are configured using googletest, and they automatically detect whether OpenMP and/or MPI are present or not and test accordingly. All the test can be run by typing "make test" from the build folder

Our report is in the file P1.7_project_report.pdf

We split the work in the following way:
1. Serial optimization: Behzad Salmassian
2. MPI implementation: Francesco Andreucci
3. OpenMP implementation: Zakaria Dahbi
4. Hybrid MPI-OpenMP: Francesco and Zakaria

The various commits begin with the name of commited. The commits up to Sunday 29/11 (included) refer to the group assignments, the ones after refer to the individual tasks.
