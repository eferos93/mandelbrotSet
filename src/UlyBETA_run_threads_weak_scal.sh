#!/bin/sh
echo "Threads version - weak scaling"
cd mandelbrotSet/src
module purge
module load intel/18.0.3.222
module load openmpi/2.1.3
mpicc -fopenmp mandelbrot_set.c -o mandelbrot_set_threads_ws.x
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

OLDIFS=$IFS; IFS=','; for procs in 2,1000 4,1414 8,2000 16,2828 24,3464 32,4000; do set -- $procs;
    export OMP_NUM_THREADS=$1
    mpirun -np 1 ./mandelbrot_set_threads_ws.x $2 $2 -2.0 -1.0 1.0 1.0 65535;
done; IFS=$OLDIFS
