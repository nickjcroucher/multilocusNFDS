#! /bin/sh

module load intel-suite
module load gsl

icc -std=c++11 -g -O2 -o freqDepSelect main.cpp functions.cpp -L /apps/gsl/1.8/lib/ -I /apps/gsl/1.8/include/  -lgsl -lgslcblas
