#!/bin/sh
cd mandelbrotSet/src
module load gnu7/7.3.1
gcc mandelbrot_set_serial.c -o mandelbrot_set_serial.x
echo "SERIAL"
./mandelbrot_set_serial.x 4000 4000 -2.0 -1.0 1.0 1.0 65535
