
echo "Hybrid version"
module load openmpi
cd mandelbrotSet/src
make

echo "PARALLEL"
OLDIFS=$IFS; IFS=','; for procs in 1,20 2,10 4,5 10,2 20,1; do IFS=","; set -- $procs;
    export OMP_NUM_THREADS=$2
    mpirun -np $1 ./mandelbrot_set.x 20000 20000 -2.0 -2.0 3.0 5.5 65535;
done;IFS=$OLDIFS
