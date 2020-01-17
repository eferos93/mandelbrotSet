cd mandelbrotSet/src
gcc mandelbrot_set_serial.c -o mandelbrot_set_serial.x -lrt
echo "SERIAL"
./mandelbrot_set_serial.x 20000 20000 -2.0 -1.0 1.0 1.0 65535
