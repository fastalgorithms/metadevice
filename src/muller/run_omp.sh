#/bin/bash

#
#  Run an OpenMP program 
#
#  Examples:
#
#  run_omp.sh 1 ./int2
#  run_omp.sh 2 ./int2
#
#  Note that the OpenMP stack size is set to 4096MB in this script.
#

PATH=.:$PATH

ulimit -s

# needed for gfortran: libgomp default stack size is too small
export OMP_STACKSIZE=4096M

(export OMP_NUM_THREADS=$1; shift; /usr/bin/time $*)

