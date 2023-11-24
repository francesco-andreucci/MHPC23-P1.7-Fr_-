#!/bin/sh
echo "BUILDING OBJECTS"
cmake -S . -B $1

echo " "
echo " "
echo "BUILDING EXECUTABLES FROM OBJECTS"
cmake --build $1