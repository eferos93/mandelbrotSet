
mandelbrot_set.x: mandelbrot_set.c mandelbrot_set_serial.x
	mpicc -std=c11 -fopenmp $< -o $@

mandelbrot_set_serial.x:  mandelbrot_set_serial.c
	gcc $< -o $@

clean:
	rm -f mandelbrot_set.x mandelbrot_set_serial.x
