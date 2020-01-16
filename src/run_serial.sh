module load openmpi
cd mandelbrotSet/src
gcc mandelbrot_set_serial.c -o mandelbrot_set_serial.x -lrt
echo "SERIAL"
./mandelbrot_set_serial.x 20000 20000 -2.0 -2.0 3.0 5.5 65535
