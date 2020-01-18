!/usr/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l walltime=01:00:00
#PBS -q devel
cd $PBS_O_WORKDIR
module purge
module load gnu
gcc mandelbrot_set_serial.c -o mandelbrot_set_serial.x -lrt
echo "SERIAL"
for input in 500 707 1000 1225 1421 1579; do
    ./mandelbrot_set_serial.x ${input} ${input} -2.0 -1.0 1.0 1.0 65535
done
