
mandelbrot_set.x: mandelbrot_set.c mandelbrot_set_serial.x
	mpicc -std=c11 -fopenmp $< -o $@

mandelbrot_set_serial.x:  mandelbrot_set_serial.c
	gcc $< -o $@ -lrt

clean:
	rm -f *.x *.o
