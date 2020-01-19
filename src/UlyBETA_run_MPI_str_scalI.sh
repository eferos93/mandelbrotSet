#!/bin/sh
echo "MPI version"
module purge
module load intel/18.0.3.222
module load openmpi/2.1.3
cd mandelbrotSet/src
mpicc mandelbrot_set.c -o mandelbrot_set_MPI.x

for procs in {2..20}; do
    mpirun -np ${procs} ./mandelbrot_set_MPI.x 3000 3000 -2.0 -1.0 1.0 1.0 65535
done
    
