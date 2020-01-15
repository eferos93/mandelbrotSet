module load openmpi
make
cd mandelbrotSet/src
echo "SERIAL"
\.mandelbrot_set_serial.x 20000 20000 -2.0 -2.0 3.0 5.5 65535