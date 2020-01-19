
#!/bin/sh
cd mandelbrotSet/src
echo "Threads version"
module purge
module load intel/18.0.3.222
module load openmpi/2.1.3
mpicc -fopenmp mandelbrot_set.c -o mandelbrot_set_threads.x

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

for procs in {2..20}; do
    export OMP_NUM_THREADS=${procs}
    mpirun -np 1 ./mandelbrot_set_threads.x 3000 3000 -2.0 -1.0 1.0 1.0 65535;
done
    