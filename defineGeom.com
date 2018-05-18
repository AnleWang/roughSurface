rm defineGeom.x
g++ -O2 defineGeom.cpp -std=c++11 -lfftw3 -o defineGeom.x >& defineGeom.e
more defineGeom.e
