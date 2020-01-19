#!/bin/sh
echo "Hybrid version"
module purge
module load intel/18.0.3.222
module load openmpi/2.1.3

cd mandelbrotSet/src
mpicc -fopenmp mandelbrot_set.c -o mandelbrot_set_hybrid.x
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=off
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

OLDIFS=$IFS; IFS=','; for procs in 1,20 4,5 10,2 20,1; do set -- $procs;
    export OMP_NUM_THREADS=$2
    mpirun -np $1 ./mandelbrot_set_hybrid.x 3000 3000 -2.0 -1.0 1.0 1.0 65535;
done; IFS=$OLDIFS
