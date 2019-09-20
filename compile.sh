#! /bin/sh
g++ -std=c++14 -pedantic-errors -Wextra  -O2 -g -o freqDepSelect  main.cpp  functions.cpp  -I /usr/local/include/ -L /usr/local/lib/ -lgsl -lgslcblas
