#!/bin/bash

cd source/qpOASES
mkdir build
cd build
cmake -G "Unix Makefiles" ..
cmake --build .

cd ..
rm -r examples
rm -r interfaces
rm -r testing
cd build
rm -r CMakeFiles
cd ..
rm Makefile 'make.mk' make_cygwin.mk make_linux.mk make_osx.mk make_windows.mk
cd ..
cd ..

make

./runner samplesMPC_N4.txt 1 1 1 10

