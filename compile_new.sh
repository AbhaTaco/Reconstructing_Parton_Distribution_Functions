#!/bin/bash

main=$1

g++ -Wall $main ./pseudo_pdf_lattice.cpp -lfftw3f -lm -o go

#./class/test_class.cpp
