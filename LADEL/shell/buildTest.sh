#!/bin/bash

#Assume this is called in the LADEL directory
curdir=`pwd`
rm -r build
#Build direcetories
if [ ! -d "build" ]; then
  mkdir build
fi

if [ ! -d "build/debug" ]; then
  mkdir build/debug
fi

if [ ! -d "build/lib" ]; then
  mkdir build/lib
fi

builddir=$curdir/build/debug
cd $builddir

cmake $curdir -DCMAKE_BUILD_TYPE=debug -DUNITTESTS=ON -DSIMPLE_COL_COUNTS=OFF -DLONG=ON -DAMD=ON
make
ctest -VV
# valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose bin/ladel_run_all_tests