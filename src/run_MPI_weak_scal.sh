echo "MPI version"
module purge
module load openmpi/1.8.3/intel/14.0
module load impi-trial/5.0.1.035
cd mandelbrotSet/src
mpiicc mandelbrot_set.c -o mandelbrot_set_MPI_ws.x
export OMP_NUM_THREADS=1

OLDIFS=$IFS; IFS=','; for procs in 2,2000 4,4000 8,8000 12,12000 16,16000 20,20000; do set -- $procs;
    mpiexec.hydra -n $1 -ppn 1 ./mandelbrot_set_MPI_ws.x $2 $2 -2.0 -2.0 3.0 5.5 65535;
done; IFS=$OLDIFS
    
