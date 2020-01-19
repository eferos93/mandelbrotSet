#!/bin/sh
echo "MPI version - weak scaling"
module purge
module load intel/18.0.3.222
module load openmpi/2.1.3
cd mandelbrotSet/src
mpicc mandelbrot_set.c -o mandelbrot_set_MPI_ws.x
export OMP_NUM_THREADS=1

OLDIFS=$IFS; IFS=','; for procs in 2,500 4,707 8,1000 12,1225 16,1421 20,1579; do set -- $procs;
    mpirun -np $1 ./mandelbrot_set_MPI_ws.x $2 $2 -2.0 -1.0 1.0 1.0 65535;
done; IFS=$OLDIFS
    
