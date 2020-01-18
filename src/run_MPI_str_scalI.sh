echo "MPI version"
module purge
module load openmpi/1.8.3/intel/14.0
module load impi-trial/5.0.1.035

cd $PBS_O_WORKDIR
mpiicc mandelbrot_set.c -o mandelbrot_set_MPI.x

for procs in {2..20}; do
    mpiexec.hydra -n ${procs} ./mandelbrot_set_MPI.x 6000 6000 -2.0 -1.0 1.0 1.0 65535
done
    
