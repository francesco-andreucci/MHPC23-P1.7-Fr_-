# CMake for the project

# version 3.16 is bundled with ubuntu20.04
cmake_minimum_required(VERSION 3.16)

# set the project name, version, and primary programming language(s) 
project(ljmd-project VERSION 0.1 LANGUAGES C CXX)

# define minimum c++ version to zakaria's mac
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(USE_MPI "Building with MPI" ON)
if(USE_MPI)
  find_package(MPI REQUIRED)
  add_compile_definitions(LJMD_MPI=yes)
  include_directories(${MPI_C_INCLUDE_PATH})
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif()



# prefetch and link openmp when found
# find_package(OpenMP COMPONENTS C REQUIRED)
# add_definitions(-DOMP_H -DOMP_V2)

# get ride of clang warnings on zakaria's mac
if(APPLE)
  add_definitions(-Qunused-arguments)
endif()

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()


#library for the c functions
add_library(mdlib  src/utilities.c src/constants.c src/input.c src/output.c src/verlet.c src/memory.c  src/force_comp.c )
if(USE_MPI)
  target_link_libraries(mdlib PUBLIC MPI::MPI_C)
endif()
target_compile_options(mdlib PRIVATE -O3 -ffast-math -fexpensive-optimizations -msse3)
# target_link_libraries(mdlib m OpenMP::OpenMP_C)
target_include_directories(mdlib PUBLIC include/)
# target_compile_options(mdlib PRIVATE "${OpenMP_CXX_FLAGS}")

# add the executable
add_executable(main src/main.cpp)
target_compile_options(main PRIVATE -O3 -ffast-math -fexpensive-optimizations -msse3)
if(USE_MPI)
  target_link_libraries(main PUBLIC mdlib m MPI::MPI_CXX)
else()
  target_link_libraries(main PUBLIC mdlib m)
endif()


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
      COMMENT "Preparing ${FILENAME} needed for input_test" VERBATIM
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
