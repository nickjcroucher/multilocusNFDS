#! /bin/sh
g++ -Werror -Wall -g -O2 -o freqDepSelect  main.cpp  functions.cpp  -I /usr/local/include/ -L /usr/local/lib/ -lgsl -lgslcblas
