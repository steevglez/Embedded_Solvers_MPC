#!/bin/bash

cd source/osqp
mkdir build
cd build
cmake -G "Unix Makefiles" ..
cmake --build .
rm -r CMakeFiles
cd ..
rm -r tests
cd ..
cd ..

make

export LD_LIBRARY_PATH=source/osqp/build/out
./runner samplesMPC_N4.txt 1 1 1 10
