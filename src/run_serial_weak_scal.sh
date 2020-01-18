cd $PBS_O_WORKDIR
gcc mandelbrot_set_serial.c -o mandelbrot_set_serial.x -lrt
echo "SERIAL"
for procs in 500 707 1000 1225 1421 1579; do
    ./mandelbrot_set_serial.x ${procs} ${procs} -2.0 -1.0 1.0 1.0 65535
done
