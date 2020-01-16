
echo "Threads version"
module purge
module load openmpi/1.8.3/intel/14.0
module load impi-trial/5.0.1.035
cd mandelbrotSet/src #used this for submitting job in the cluster
mpiicc -fopenmp mandelbrot_set.c -o mandelbrot_set_threads_ws.x
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=off
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

OLDIFS=$IFS; IFS=','; for procs in 2,2000 4,4000 8,8000 12,12000 16,16000 20,20000; do set -- $procs;
    export OMP_NUM_THREADS=$1
    mpiexec.hydra -n 1 -ppn 1 ./mandelbrot_set_threads_ws.x $2 $2 -2.0 -2.0 3.0 5.5 65535;
done; IFS=$OLDIFS