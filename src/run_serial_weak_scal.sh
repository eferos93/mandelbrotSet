cd $PBS_O_WORKDIR
gcc mandelbrot_set_serial.c -o mandelbrot_set_serial.x -lrt
echo "SERIAL"
for procs in 1500 2121 3000 3674 4237 4737; do
    ./mandelbrot_set_serial.x ${procs} ${procs} -2.0 -1.0 1.0 1.0 65535
done
