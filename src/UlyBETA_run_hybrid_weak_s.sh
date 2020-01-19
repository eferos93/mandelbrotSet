#!/bin/sh
echo "Hybrid version"
module purge
module load intel/18.0.3.222
module load openmpi/2.1.3

cd mandelbrotSet/src
mpicc -fopenmp mandelbrot_set.c -o mandelbrot_set_hybrid_ws.x
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

OLDIFS=$IFS; IFS=','; for procs in 1,2,1000 2,4,1414 2,8,2000 3,8,2828 5,6,3464 5,8,4000; do set -- $procs;
    export OMP_NUM_THREADS=$2
    mpirun -np $1 ./mandelbrot_set_hybrid_ws.x $3 $3 -2.0 -1.0 1.0 1.0 65535;
done; IFS=$OLDIFS
