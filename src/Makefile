#USE THIS ONLY IN LOCAL MACHINE, FOR TESTING ON THE CLUSTER, SEE THE .SH FILES
mandelbrot_set.x: mandelbrot_set.c mandelbrot_set_serial.x
	mpicc -fopenmp $< -o $@

mandelbrot_set_serial.x:  mandelbrot_set_serial.c
	gcc $< -o $@ -lrt

clean:
	rm -f *.x *.o
