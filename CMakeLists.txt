# CMake for the project

# version 3.16 is bundled with ubuntu20.04
cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE STRING "C++ compiler" FORCE)
set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE STRING "C compiler" FORCE)
# set the project name, version, and primary programming language(s) 
project(ljmd-project VERSION 0.1 LANGUAGES C CXX)


# define minimum c++ version to zakaria's mac
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# get ride of clang warnings on zakaria's mac
# if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#   add_definitions(-Qunused-arguments)
# endif()

set(LJMD_CUSTOM_COMPILER_FLAGS "-O3 -Wall -g -Ofast")


if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(LJMD_CUSTOM_COMPILER_FLAGS "${LJMD_CUSTOM_COMPILER_FLAGS} -Ofast -Xclang -fopenmp")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(LJMD_CUSTOM_COMPILER_FLAGS "${LJMD_CUSTOM_COMPILER_FLAGS} -fexpensive-optimizations")
endif()

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()

# turn ON/OFF OMP support
option(USE_OMP "Building with OPENMP" ON)
if(USE_OMP)
  find_package(OpenMP REQUIRED)
  #add_definitions(-fopenmp)
  add_compile_definitions(LJMD_OPENMP=yes)
endif()

# turn ON/OFF MPI support
option(USE_MPI "Building with MPI" ON)
if(USE_MPI)
  find_package(MPI REQUIRED)
  add_compile_definitions(LJMD_MPI=yes)
  include_directories(${MPI_C_INCLUDE_PATH})
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif()

#library for the c functions
add_library(mdlib  src/utilities.c src/constants.c src/input.c src/output.c src/verlet.c src/memory.c  src/force_comp.c )

# compile mdlib with OMP support if used
if(USE_OMP)
  target_link_libraries(mdlib OpenMP::OpenMP_C)
endif()

# compile mdlib with MPI support if used
# if(USE_MPI)
#   target_link_libraries(mdlib PUBLIC MPI::MPI_C)
# endif()

# compiling mdlib with optimization flags
target_compile_options(mdlib PRIVATE -O3 -g -ffast-math -fexpensive-optimizations -msse3 -Wall)
target_include_directories(mdlib PUBLIC include/)

# add the executable
add_executable(main src/main.cpp)

target_compile_options(main PRIVATE -O3 -ffast-math -fexpensive-optimizations -msse3 -Wall)
if(USE_MPI)
  target_link_libraries(main PUBLIC mdlib m MPI::MPI_CXX)
else()
  target_link_libraries(main PUBLIC mdlib m)
endif()

 if(USE_OMP)
   target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)

   if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
   target_compile_options(main PUBLIC -Xclang -fopenmp)
   else()
    target_compile_options(main PUBLIC -fopenmp)
   endif()

 endif()


# copy input files from examples to build whenever created
set(INPUT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/examples)
set(WILDCARD_PATTERN_INPUT "*")
file(GLOB FILES_TO_COPY ${INPUT_FOLDER}/${WILDCARD_PATTERN_INPUT})
list(FILTER FILES_TO_COPY EXCLUDE REGEX ".*Makefile$")
foreach(FILE ${FILES_TO_COPY})
get_filename_component(FILENAME ${FILE} NAME)
add_custom_command(TARGET main POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${FILE} $<TARGET_FILE_DIR:main>/${FILENAME}
    COMMENT "Copied input file: ${FILENAME} to build -> needed for I/O testing" VERBATIM
)
endforeach()

# To enable of disable testing
option(ENABLE_TESTING "Enable building unit tests" ON)
if(ENABLE_TESTING)
  enable_testing()
  enable_language(CXX)
  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  include(GTest)

  # gtest for test_force
  add_executable(test_force tests/test_forces.cpp)
  target_include_directories(test_force PRIVATE include/)

  if(USE_MPI)
  target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest mdlib MPI::MPI_CXX)
  else()
  target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest mdlib)
  endif()

  if(USE_OMP)
  target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest mdlib OpenMP::OpenMP_CXX)
  
  if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  target_compile_options(test_force PUBLIC -Xclang -fopenmp)
  else()
   target_compile_options(test_force PUBLIC -fopenmp)
  endif()
  
  # target_compile_options(test_force PUBLIC  -fopenmp)
  else()
  target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest mdlib)
  endif()

  # if(USE_OMP)
  # target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest mdlib  OpenMP::OpenMP_CXX)
  # target_compile_options(test_force PUBLIC -Xclang -fopenmp)
  # endif()


  add_test(NAME Force COMMAND test_force)
# gtest for test_ekin
  add_executable(test_ekin tests/test_ekin.cpp)
  target_include_directories(test_ekin PRIVATE include/)
  target_link_libraries(test_ekin PRIVATE GTest::GTestMain GTest::GTest mdlib m)
  add_test(NAME EkinTest COMMAND test_ekin)
# gtest for test_Verlet
  add_executable(test_verlet tests/test_verlet.cpp)
  target_include_directories(test_verlet PRIVATE include/)
  target_link_libraries(test_verlet PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME VerletTest COMMAND test_verlet)
  # gtest for test_input
  # copy ekin_test for input_test
  add_executable(test_input tests/test_input.cpp)
  set(INPUT_TEST_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/tests/input_test)
  set(WILDCARD_PATTERN "ekin_test.*")
  file(GLOB FILES_TO_COPY ${INPUT_TEST_FOLDER}/${WILDCARD_PATTERN})
  foreach(FILE ${FILES_TO_COPY})
  get_filename_component(FILENAME ${FILE} NAME)
  add_custom_command(TARGET test_input POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${FILE} $<TARGET_FILE_DIR:test_input>/${FILENAME}
      COMMENT "Copied ${FILENAME} to build -> needed for input_test" VERBATIM
  )
endforeach()
# run gtest for input_test
  target_include_directories(test_input PRIVATE include/)
  if(USE_MPI)
  target_link_libraries(test_input PRIVATE GTest::GTestMain GTest::GTest mdlib MPI::MPI_CXX)
  else()
  target_link_libraries(test_input PRIVATE GTest::GTestMain GTest::GTest mdlib)
  endif()
  add_test(NAME InputTest
  COMMAND sh -c "exec $<TARGET_FILE:test_input> < ekin_test.inp")
endif()
