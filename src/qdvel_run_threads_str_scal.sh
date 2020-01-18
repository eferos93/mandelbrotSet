
!/usr/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l walltime=01:00:00
#PBS -q devel
cd $PBS_O_WORKDIR
echo "Threads version"
module purge
module load openmpi/1.8.3/intel/14.0
module load impi-trial/5.0.1.035
mpiicc -fopenmp mandelbrot_set.c -o mandelbrot_set_threads.x
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=off
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

for procs in {2..20}; do
    export OMP_NUM_THREADS=${procs}
    mpiexec.hydra -n 1 ./mandelbrot_set_threads.x 2000 2000 -2.0 -1.0 1.0 1.0 65535;
done
    