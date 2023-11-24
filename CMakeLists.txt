# CMake for the project

# version 3.16 is bundled with ubuntu20.04
cmake_minimum_required(VERSION 3.16)

# set the project name, version, and primary programming language(s) 
project(ljmd-project VERSION 0.1 LANGUAGES C CXX)

#library for the c functions
add_library(mdlib  src/utilities.c src/constants.c src/input.c src/output.c src/verlet.c src/cleanup.c  src/force_comp.c )
target_link_libraries(mdlib  m)
target_include_directories(mdlib PUBLIC include/)




# add the executable
add_executable(main src/main.cpp)
target_link_libraries(main mdlib m)
