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
#include <mpi.h>
#include <stdbool.h>

//mpi tags
#define DATA 0
#define RESULT 1
#define TERMINATE 2
#define REMAINDER 3
#define MASTER 0
#define LOAD_BALANCE_TIME

unsigned int I_max, N_x, N_y, job_width, job_remainder, job_remainder_size, job_size, *image, *result;
int pid, world_size;
double x_L, x_R, y_L, y_R;
double d_x, d_y;

struct complex
{
    double real;
    double imag;
};

/**
 * Function that given a a complex number c and and integer,
 * computes if c belongs to the Mandelbrot Set and returns
 * the counter used in the loop 
 */
int cal_pixel(struct complex c, unsigned int max_iter) {
 int count=0;
 struct complex z = c;
 double temp;

 while ((z.real * z.real + z.imag * z.imag <= 4.0) && (count < max_iter)) {
    temp = z.real * z.real - z.imag * z.imag + c.real;
    z.imag = 2 * z.real * z.imag + c.imag;
    z.real = temp;
    count++;
 }

 return count;
}


/**
 * prints a grayscale image
 */
void write_pgm_image( void *image, int maxval, unsigned int xsize, unsigned int ysize, const char *image_name)
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
 * methods representing the work to be carried out by a slave
 * start represent the starting index from which the elements computed
 * will be placed in the final array of pixels
 * result is a pointer to an array of fixed length which first element
 * will contain the before mentioned start and the elements computed
 */
void _worker(unsigned int start, unsigned int* result, unsigned int work_amount)
{
    struct complex c;
    result[0] = start;

    #pragma omp parallel for schedule(dynamic, 10) private(c) collapse(2)
    for(int j=0; j < N_y; j++) {
        for(int i = 0; i < work_amount; i++) {
            c.real = (i + start) * d_x + x_L;
            c.imag = j * d_y + y_L;
            result[j * work_amount + i + 1] = cal_pixel(c, I_max);
        }
    }
}

/**
 * method that will be executed by the master only, which will act
 * as a manager that redistributes the work
 */
void master()
{
    if (world_size == 1) {
        _worker(MASTER, result);
        return;
    }

    bool remainder = false;
    unsigned int iterat, start;
    double start_wtime = MPI_Wtime();
    MPI_Status stat;
    unsigned int actives = 1, jobs = 0;

    for (; actives < world_size && jobs < N_x; actives++, jobs += job_width) {
        int tag;
        if(jobs + job_width > N_x) tag = REMAINDER; else tag = DATA;
        
        MPI_Send(jobs, 1, MPI_UNSIGNED, actives, tag, MPI_COMM_WORLD);
    }
    
    do {
        result = (unsigned int*) malloc((job_size + 1) * sizeof(int));
        MPI_Recv(result, job_size+1, MPI_UNSIGNED, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
        if(stat.MPI_TAG == DATA) {
            iterat = job_size;
        } else {
            iterat = job_remainder_size;
        }
        
        
        int slave = stat.MPI_SOURCE;
        start = result[0];
        //no need to parallelize this loop, just 20 iterations
        for(int i = 1; start < start + iterat; i++, start++) {
            image[start] = result[i];
        }
        
        actives--;
        if (jobs < N_x) { 
            if (jobs + job_width > N_x) {
                //if we enter in this branch, we are issuing the last job
                //and N_y is not divisible by job_width (i.e.  0 < job_remainder < 20)
                //therefore we send a message with a different tag
                printf("1\n");
                MPI_Send(jobs, 1, MPI_UNSIGNED, slave, REMAINDER, MPI_COMM_WORLD);
                jobs += job_remainder;
            } else {
                printf("2\n");
                MPI_Send(jobs, 1, MPI_UNSIGNED, slave, DATA, MPI_COMM_WORLD);
                jobs += job_width;
            }

            actives++;
        } else { 
            MPI_Send(jobs, 1, MPI_UNSIGNED, slave, TERMINATE, MPI_COMM_WORLD);
        }
    } while (actives > 1);

    printf("Total Wtime: %f", MPI_Wtime() - start_wtime);
    write_pgm_image((void*) image, I_max, N_x, N_y, "mandelbrot_set.pgm");
}

/**
 * method used by slave processes, that will carry out the work and 
 * reply back to the master with the results until the the 
 * MPI_TAG == DATA
 */
void slave()
{
    MPI_Status stat;
    unsigned int col;
    MPI_Recv(col, 1, MPI_UNSIGNED, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    printf("Slave %d col: %d received", pid, col);
    result = (unsigned int*) malloc((job_size + 1) * sizeof(int));

    while (stat.MPI_TAG != TERMINATE) 
    {
        int tag;
        if(stat.MPI_TAG == DATA) {
            _worker(col, result, job_width);
            tag = DATA;
        } else {
            //TAG == REMAINDER
            _worker(col, result, job_remainder);
            tag = REMAINDER;
        }
        MPI_Send(result, job_size+1, MPI_UNSIGNED, MASTER, tag, MPI_COMM_WORLD);
        MPI_Recv(col, 1, MPI_UNSIGNED, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    }
}

/**
 * Initialise enviroment
 */
void initial_env(int argc, char** argv) {
    N_x = atoi(argv[1]), N_y = atoi(argv[2]);
    x_L = atof(argv[3]), y_L = atof(argv[4]);
    x_R = atof(argv[5]), y_R = atof(argv[6]);
    I_max = atoi(argv[7]);

    if(N_y < 250 || N_y < 250) {
        printf("N_y and N_x must be at least 250");
        exit(0);
    }
    d_x = (x_R - x_L) / N_x;
    d_y = (y_R - y_L) / N_y;
}

/**
 * Initialise MPI environment
 */
void initial_MPI_env(int argc, char** argv)
{
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    if(pid == MASTER) {
        image = (unsigned int*) malloc(N_y * N_x * sizeof(int));
    }
    
    job_width = world_size == 1 ? N_x : 20;
    if (world_size > 1) {
        job_remainder = (N_x % 20);
        job_remainder_size = job_remainder * N_y;
    }
    job_size = job_width * N_y;
}

void start()
{
    pid == MASTER ? master() : slave();
}

int main(int argc, char** argv) {
    initial_env(argc, argv);
    printf("Init done\n");
    initial_MPI_env(argc, argv);
    printf("MPI init done\n");
    //this matrix is 20x20, no need to allocate it in the heap
    //int thread_time[world_size][omp_get_num_threads];
    start();
    MPI_Finalize();
    return 0;
}