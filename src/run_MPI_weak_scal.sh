echo "MPI version - weak scaling"
module purge
module load openmpi/1.8.3/intel/14.0
module load impi-trial/5.0.1.035
cd mandelbrotSet/src
mpiicc mandelbrot_set.c -o mandelbrot_set_MPI_ws.x
export OMP_NUM_THREADS=1

OLDIFS=$IFS; IFS=','; for procs in 2,1500 4,2121 8,3000 12,3674 16,4237 20,4737; do set -- $procs;
    mpiexec.hydra -n $1 -ppn 1 ./mandelbrot_set_MPI_ws.x $2 $2 -2.0 -1.0 1.0 1.0 65535;
done; IFS=$OLDIFS
    
