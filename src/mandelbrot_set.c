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
 * Function that given a complex number c and and integer,
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
    printf("Slave %d in _worker; start: %u; work_amount: %u\n", pid, start, work_amount);
    #pragma omp parallel for schedule(dynamic, 10) private(c) collapse(2)
    for(int i = 0; i < work_amount; i++) {
        for(int j = 0; j < N_y; j++) {
            c.real = (i + start) * d_x + x_L;
            c.imag = j * d_y + y_L;
            result[j * work_amount + i + 1] = compute_pixel(c, I_max);
        }
    }
}

/**
 * method that will be executed by the master only, which will act
 * as a manager that redistributes the work
 */
void master()
{
    printf("Master started\n");
    if (world_size == 1) {
        _worker(MASTER, result, N_x);
        return;
    }
    unsigned int iterat, start;
    MPI_Status stat;
    unsigned int actives = 1, jobs = 0;
    int tag = DATA;

    double start_wtime = MPI_Wtime();

    for (; actives < world_size && jobs < N_x; actives++, jobs += job_width) {
        //if N_x % 20 != 0 then, for the last message to be sent, we use a different tag
        if(jobs + job_width > N_x) tag = REMAINDER;
        printf("Sending %u offset to %u slave with tag %d\n", jobs, actives, tag);
        MPI_Send(&jobs, 1, MPI_UNSIGNED, actives, tag, MPI_COMM_WORLD);
    }
    printf("jobs distributed\n");
    do {
        
        result = (unsigned int*) malloc((job_size + 1) * sizeof(unsigned int));
        
        if (result == NULL) 
        {
            printf("Not able to allocate result on MASTER");
            exit(1);
        }

        MPI_Recv(result, job_size+1, MPI_UNSIGNED, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
        
        iterat = stat.MPI_TAG == DATA ? job_size : job_remainder_size;
        
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
                MPI_Send(&jobs, 1, MPI_UNSIGNED, slave, REMAINDER, MPI_COMM_WORLD);
                jobs += job_remainder;
            } else {
                printf("2\n");
                MPI_Send(&jobs, 1, MPI_UNSIGNED, slave, DATA, MPI_COMM_WORLD);
                jobs += job_width;
            }

            actives++;
        } else { 
            MPI_Send(&jobs, 1, MPI_UNSIGNED, slave, TERMINATE, MPI_COMM_WORLD);
        }
        
        free(result);
    } while (actives > 1);

    printf("Total Wtime: %f\n", MPI_Wtime() - start_wtime);
    write_pgm_image(image, I_max, N_x, N_y, "mandelbrot_set_parallel.pgm");
}

/**
 * method used by slave processes, that will carry out the work and 
 * reply back to the master with the results until the the 
 * MPI_TAG == DATA
 */
void slave()
{
    printf("Slave %d started\n",pid);
    MPI_Status stat;
    unsigned int column_offset;
    int tag;
    MPI_Recv(&column_offset, 1, MPI_UNSIGNED, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    tag = stat.MPI_TAG;
    printf("Slave %d column_offset: %d received\n", pid, column_offset);
    result = (unsigned int*) malloc((job_size + 1) * sizeof(unsigned int));
    if (result == NULL)
    {
        printf("Malloc error in slave %d\n", pid);
    }
    
    while (tag != TERMINATE) 
    {
        
        if(tag == DATA) {
            printf("Slave: %d; issung _worker with col_off: %u; job_width: %u\n", pid, column_offset, job_width);
            _worker(column_offset, result, job_width);
        } else {
            //TAG == REMAINDER
            printf("Slave %d; issung _worker with col_off: %u; job_rem: %u\n", pid, column_offset, job_remainder);
            _worker(column_offset, result, job_remainder);
        }
        MPI_Send(result, job_size+1, MPI_UNSIGNED, MASTER, tag, MPI_COMM_WORLD);
        MPI_Recv(&column_offset, 1, MPI_UNSIGNED, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
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

    if(N_y < 250 || N_y < 250) {
        printf("N_y and N_x must be at least 250");
        exit(0);
    }
    if(N_y != N_x) {
        printf("N_x and N_x must be equal");
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
        image = (unsigned int*) malloc(N_y * N_x * sizeof(unsigned int));
    }
    
    job_width = world_size == 1 ? N_x : 20;
    if (world_size > 1) {
        job_remainder = (N_x % 20);
        job_remainder_size = job_remainder * N_y;
    }
    job_size = job_width * N_y;
    printf("job_remainder: %u; job_width: %u\n", job_remainder, job_width);
}

int main(int argc, char** argv) {
    initial_env(argc, argv);
    printf("Init done\n");
    initial_MPI_env(argc, argv);
    printf("pid %d: MPI init done\n", pid);
    pid == MASTER ? master() : slave();
    
    MPI_Finalize();
    return 0;
}
