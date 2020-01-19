#!/bin/sh
cd mandelbrotSet/src
module load gnu7/7.3.1
gcc mandelbrot_set_serial.c -o mandelbrot_set_serial.x
echo "SERIAL"
for input in 500 707 1000 1225 1421 1579; do
    ./mandelbrot_set_serial.x ${input} ${input} -2.0 -1.0 1.0 1.0 65535
done
