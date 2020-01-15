echo "MPI version"
module load openmpi
cd mandelbrotSet/src
mpicc -std=c11 mandelbrot_set.c -o mandelbrot_set.x

ECHO "PARALLEL"
for procs in {2..20}; do
    mpirun -np ${procs} ./mandelbrot_set.x 20000 20000 -2.0 -2.0 3.0 5.5 65535
done
    
