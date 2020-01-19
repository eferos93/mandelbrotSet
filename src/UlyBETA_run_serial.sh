#!/bin/sh
cd mandelbrotSet/src
gcc mandelbrot_set_serial.c -o mandelbrot_set_serial.x
echo "SERIAL"
./mandelbrot_set_serial.x 2000 2000 -2.0 -1.0 1.0 1.0 65535
