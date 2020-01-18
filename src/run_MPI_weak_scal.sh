echo "MPI version - weak scaling"
module purge
module load openmpi/1.8.3/intel/14.0
module load impi-trial/5.0.1.035
cd $PBS_O_WORKDIR
mpiicc mandelbrot_set.c -o mandelbrot_set_MPI_ws.x
export OMP_NUM_THREADS=1

OLDIFS=$IFS; IFS=','; for procs in 2,500 4,707 8,1000 12,1225 16,1421 20,1579; do set -- $procs;
    mpiexec.hydra -n $1 -ppn 1 ./mandelbrot_set_MPI_ws.x $2 $2 -2.0 -1.0 1.0 1.0 65535;
done; IFS=$OLDIFS
    
