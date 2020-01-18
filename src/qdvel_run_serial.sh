!/usr/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l walltime=01:00:00
#PBS -q devel
cd $PBS_O_WORKDIR
module purge
module load gnu
gcc mandelbrot_set_serial.c -o mandelbrot_set_serial.x -lrt
echo "SERIAL"
./mandelbrot_set_serial.x 2000 2000 -2.0 -1.0 1.0 1.0 65535
