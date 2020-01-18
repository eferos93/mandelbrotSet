!/usr/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l walltime=01:00:00
#PBS -q devel
cd $PBS_O_WORKDIR
echo "MPI version"
module purge
module load openmpi/1.8.3/intel/14.0
module load impi-trial/5.0.1.035
mpiicc mandelbrot_set.c -o mandelbrot_set_MPI.x

for procs in {2..20}; do
    mpiexec.hydra -n ${procs} ./mandelbrot_set_MPI.x 2000 2000 -2.0 -1.0 1.0 1.0 65535
done
    
