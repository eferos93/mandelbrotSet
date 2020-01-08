
mandelbrot_set.x: mandelbrotSet.c
	mpicc -fopenmp <$ -o $@

mandelbrot_set_serial.x: mandelbrot_set_serial.c
	gcc <$ -o $@