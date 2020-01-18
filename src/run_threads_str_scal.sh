
echo "Threads version"
module purge
module load openmpi/1.8.3/intel/14.0
module load impi-trial/5.0.1.035
cd mandelbrotSet/src #used this for submitting job in the cluster
mpiicc -fopenmp mandelbrot_set.c -o mandelbrot_set_threads.x
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=off
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

for procs in {2..20}; do
    export OMP_NUM_THREADS=${procs}
    mpiexec.hydra -n 1 ./mandelbrot_set_threads.x 6000 6000 -2.0 -1.0 1.0 1.0 65535;
done
    