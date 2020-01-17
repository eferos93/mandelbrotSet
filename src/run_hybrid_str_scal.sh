
echo "Hybrid version"
module purge
module load openmpi/1.8.3/intel/14.0
module load impi-trial/5.0.1.035
cd mandelbrotSet/src
mpiicc -fopenmp mandelbrot_set.c -o mandelbrot_set_hybrid.x
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=off
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

OLDIFS=$IFS; IFS=','; for procs in 1,20 4,6 8,3 20,1; do set -- $procs;
    export OMP_NUM_THREADS=$2
    mpiexec.hydra -n $1 -ppn 1 ./mandelbrot_set_hybrid.x 20000 20000 -2.0 -2.0 3.0 5.5 65535;
done; IFS=$OLDIFS
