cd mandelbrotSet/src
gcc mandelbrot_set_serial.c -o mandelbrot_set_serial.x -lrt
echo "SERIAL"
for procs in 4000 5657 8000 9798 11300 12633; do
    ./mandelbrot_set_serial.x ${procs} ${procs} -2.0 -1.0 1.0 1.0 65535
done
