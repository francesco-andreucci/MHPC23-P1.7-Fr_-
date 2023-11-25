# CMake for the project

# version 3.16 is bundled with ubuntu20.04
cmake_minimum_required(VERSION 3.16)

# set the project name, version, and primary programming language(s) 
project(ljmd-project VERSION 0.1 LANGUAGES C CXX)

# define minimum c++ version to zakaria's mac
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()


#library for the c functions
add_library(mdlib  src/utilities.c src/constants.c src/input.c src/output.c src/verlet.c src/memory.c  src/force_comp.c )
target_link_libraries(mdlib  m)
target_include_directories(mdlib PUBLIC include/)




# add the executable
add_executable(main src/main.cpp)
target_link_libraries(main mdlib m)


option(ENABLE_TESTING "Enable building unit tests" ON)
if(ENABLE_TESTING)
  enable_testing()

  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  include(GTest)

  # gtest for test_force
  add_executable(test_force tests/test_forces.cpp)
  target_include_directories(test_force PRIVATE include/)
  target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest mdlib)
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
endif()

