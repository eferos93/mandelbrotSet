
echo "Threads version - weak scaling"
module purge
module load openmpi/1.8.3/intel/14.0
module load impi-trial/5.0.1.035
cd mandelbrotSet/src #used this for submitting job in the cluster
mpiicc -fopenmp mandelbrot_set.c -o mandelbrot_set_threads_ws.x
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=off
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

OLDIFS=$IFS; IFS=','; for procs in 2,1500 4,2121 8,3000 12,3674 16,4237 20,4737; do set -- $procs;
    export OMP_NUM_THREADS=$1
    mpiexec.hydra -n 1 -ppn 1 ./mandelbrot_set_threads_ws.x $2 $2 -2.0 -1.0 1.0 1.0 65535;
done; IFS=$OLDIFS
