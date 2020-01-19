#!/bin/sh
module load gnu7/7.3.1
echo "SERIAL"
for input in 1000 1414 2000 2828 3464 4000; do
    ./mandelbrot_set_serial.x ${input} ${input} -2.0 -1.0 1.0 1.0 65535
done
