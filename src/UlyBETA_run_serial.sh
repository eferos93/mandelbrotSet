#!/bin/sh
module load gnu7/7.3.1
echo "SERIAL"
./mandelbrot_set_serial.x 4000 4000 -2.0 -1.0 1.0 1.0 65535
