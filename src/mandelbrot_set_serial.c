/**
 * Copyright (c) 2020 Eros Fabrici - eros.fabrici@gmail.com

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

#define _GNU_SOURCE

#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
	                  (double)ts.tv_nsec * 1e-9)

struct timespec ts;
unsigned int N_x, N_y, I_max, *image;
double x_L, x_R, y_L, y_R;
double d_x, d_y;

struct complex
{
    double real;
    double imag;
};

void write_pgm_image(void *image, unsigned int maxval,
                     unsigned int xsize, unsigned int ysize, const char *image_name);

unsigned int compute_pixel(struct complex c, unsigned int max_iter);
void _worker(unsigned int* result);
void initial_env(int argc, char** argv);


int main(int argc, char** argv) {
    initial_env(argc, argv);
    image = (unsigned int*) malloc(N_x * N_y * sizeof(unsigned int));
    double start = CPU_TIME;
    _worker(image);
    printf("Total time: %f", CPU_TIME - start);
    write_pgm_image(image, I_max, N_x, N_y, "mandelbrot_set.pgm");
    return 0;
}

/**
 * prints a grayscale image
 */
void write_pgm_image(void *image, unsigned int maxval,
                     unsigned int xsize, unsigned int ysize, const char *image_name)
{
  FILE* image_file; 
  image_file = fopen(image_name, "w"); 
  
  // Writing header
  // The header's format is as follows, all in ASCII.
  // "whitespace" is either a blank or a TAB or a CF or a LF
  // - The Magic Number (see below the magic numbers)
  // - the image's width
  // - the height
  // - a white space
  // - the image's height
  // - a whitespace
  // - the maximum color value, which must be between 0 and 65535
  //
  // if he maximum color value is in the range [0-255], then
  // a pixel will be expressed by a single byte; if the maximum is
  // larger than 255, then 2 bytes will be needed for each pixel
  //

  int color_depth = 1+((maxval>>8)>0);       // 1 if maxval < 256, 2 otherwise

  fprintf(image_file, "P5\n%d %d\n%d\n", xsize, ysize, maxval);
  
  // Writing file
  fwrite( image, color_depth, xsize*ysize, image_file);  

  fclose(image_file); 
  return ;

  /* ---------------------------------------------------------------

     TYPE    MAGIC NUM     EXTENSION   COLOR RANGE
           ASCII  BINARY

     PBM   P1     P4       .pbm        [0-1]
     PGM   P2     P5       .pgm        [0-255]
     PPM   P3     P6       .ppm        [0-2^16[
  
  ------------------------------------------------------------------ */
}

/**
 * Function that given a a complex number c and and integer,
 * computes if c belongs to the Mandelbrot Set and returns
 * the counter used in the loop 
 */
unsigned int compute_pixel(struct complex c, unsigned int max_iter) 
{
    unsigned int count=0;
    struct complex z;
    z.real = 0.0;
    z.imag = 0.0;
    double temp;

    while ((z.real * z.real + z.imag * z.imag < 4.0) && (count < max_iter)) {
        temp = z.real * z.real - z.imag * z.imag + c.real;
        z.imag = 2.0 * z.real * z.imag + c.imag;
        z.real = temp;
        count++;
    }

    return count;
}


void _worker(unsigned int* buffer)
{
    struct complex c;

    for (int i = 0; i < N_y; i++)
    {
        c.imag = y_L + i * d_y;
        for (int j = 0; j < N_x; j++)
        {
            c.real = x_L + j * d_x;
            *(buffer + (i * N_x + j)) = compute_pixel(c, I_max);
        }
        
    }
    
}

/**
 * Initialise enviroment
 */
void initial_env(int argc, char** argv) 
{
    N_x = atoi(argv[1]), N_y = atoi(argv[2]);
    x_L = atof(argv[3]), y_L = atof(argv[4]);
    x_R = atof(argv[5]), y_R = atof(argv[6]);
    I_max = atoi(argv[7]);
    
    d_x = (x_R - x_L) / N_x;
    d_y = (y_R - y_L) / N_y;
    printf("N_x: %d; N_y: %d; x_L: %f; y_L: %f; x_R: %f; y_R: %f; I_max: %d; d_x: %f; d_y: %f;",
            N_x, N_y, x_L, y_L, x_R, y_R, I_max, d_x, d_y);
}