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
#include <string.h>

//mpi tags
#define DATA 0
#define RESULT 1
#define TERMINATE 2
#define REMAINDER 3
#define MASTER 0
#define LOAD_BALANCE_TIME

unsigned int I_max, N_x, N_y, job_height, job_remainder, job_remainder_size, job_size, header_size;
int pid, world_size;
double x_L, x_R, y_L, y_R;
double d_x, d_y;
char* header;
MPI_File output_file;
MPI_Status file_stat;

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
 * methods representing the work to be carried out by a slave
 * start represent the starting index from which the elements computed
 * will be placed in the final array of pixels
 * result is a pointer to an array of fixed length which first element
 * will contain the before mentioned start and the elements computed
 */
void _worker(unsigned int start, unsigned int work_amount)
{
    struct complex c;
    unsigned int* buffer = (unsigned int*) malloc(work_amount * N_x * sizeof(unsigned int));
    printf("Slave %d in _worker; start: %u; work_amount: %u\n", pid, start, work_amount);
    //#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 10) private(c) collapse(2)
    //#endif
    for(unsigned int i = 0; i < work_amount; i++) {
        for(unsigned int j = 0; j < N_x; j++) {
            c.imag = y_L + (i + start) * d_y;
            c.real = x_L + j * d_x;
            *(buffer + (i * N_x + j)) = compute_pixel(c, I_max);
            //printf("Slave: %d; result[%u] = %u\n", pid, i*work_amount+j+1, *(result + (j*work_amount+i+1)));
        }
    }

    MPI_File_write_at(output_file, header_size * sizeof(char) + start * N_x * sizeof(unsigned int), buffer,
                      N_y * work_amount, MPI_UNSIGNED, &file_stat);
    free(buffer);
}

/**
 * method that will be executed by the master only, which will act
 * as a manager that redistributes the work
 */
void master()
{
    MPI_Status stat, file_stat;
    MPI_File_open(MPI_COMM_SELF, "mandelbrot_set_parallel.pgm",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
    MPI_File_write(output_file, header, header_size, MPI_CHAR, &file_stat);
    MPI_File_close(&output_file);
    MPI_File_open(MPI_COMM_WORLD, "mandelbrot_set_parallel.pgm",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
    if (world_size == 1) {
        _worker(MASTER, N_y);
        return;
    }

    unsigned int start;
    unsigned int working_slaves = 1, jobs = 0;
    int tag;

    double start_wtime = MPI_Wtime();

    for (; working_slaves < world_size && jobs < N_x; working_slaves++, jobs += job_height) {
        tag = (jobs + job_height > N_y) ? REMAINDER : DATA;
        printf("Sending %u offset to %u slave with tag %d\n", jobs, working_slaves, tag);
        MPI_Send(&jobs, 1, MPI_UNSIGNED, working_slaves, tag, MPI_COMM_WORLD);
    }

    do {
        short done;
        MPI_Recv(&done, 1, MPI_SHORT, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
        int slave = stat.MPI_SOURCE;
        working_slaves--;

        if (jobs < N_y) { 
            if (jobs + job_height > N_y) {
                //if we enter in this branch, we are issuing the last job
                //and N_x is not divisible by job_height (i.e.  0 < job_remainder < 20)
                //therefore we send a message with a different tag
                MPI_Send(&jobs, 1, MPI_UNSIGNED, slave, REMAINDER, MPI_COMM_WORLD);
                jobs += job_remainder;
            } else {
                MPI_Send(&jobs, 1, MPI_UNSIGNED, slave, DATA, MPI_COMM_WORLD);
                jobs += job_height;
            }

            working_slaves++;
        } else { 
            MPI_Send(&jobs, 1, MPI_UNSIGNED, slave, TERMINATE, MPI_COMM_WORLD);
        }
        
    } while (working_slaves > 1);

    printf("Total Wtime: %f\n", MPI_Wtime() - start_wtime);
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
    unsigned int row_offset;
    MPI_File_open(MPI_COMM_WORLD, "mandelbrot_set_parallel.pgm",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);

    MPI_Recv(&row_offset, 1, MPI_UNSIGNED, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    printf("Slave %d column_offset: %d received\n", pid, row_offset);

    while (stat.MPI_TAG != TERMINATE) 
    {
        unsigned int temp_job_height = stat.MPI_TAG == DATA ? job_height : job_remainder; 
        printf("Slave: %d; issung _worker with col_off: %u; job_width: %u\n",
                pid, row_offset, temp_job_height);
        _worker(row_offset, temp_job_height);
        printf("_worker in slave %d has finished, now sending\n", pid);
        short done = 1;
        MPI_Send(&done, 1, MPI_SHORT, MASTER, stat.MPI_TAG, MPI_COMM_WORLD);
        MPI_Recv(&row_offset, 1, MPI_UNSIGNED, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
        if (stat.MPI_TAG == TERMINATE) 
        {
            printf("Slave %d received a message with tag = TERMINATE\n", pid);
        }
    }
}

/**
 * Initialise enviroment
 */
void init_env(int argc, char** argv) 
{
    N_x = atoi(argv[1]), N_y = atoi(argv[2]);
    x_L = atof(argv[3]), y_L = atof(argv[4]);
    x_R = atof(argv[5]), y_R = atof(argv[6]);
    I_max = atoi(argv[7]);
    
    d_x = (x_R - x_L) / N_x;
    d_y = (y_R - y_L) / N_y;
}

/**
 * Initialise MPI environment
 */
void init_MPI_env(int argc, char** argv)
{
    int mpi_provided_thread_level;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,
                     &mpi_provided_thread_level);
    if (mpi_provided_thread_level < MPI_THREAD_FUNNELED) 
    {
        printf("a problem arise when asking for MPI_THREAD_FUNNELED level\n");
        MPI_Finalize();
        exit( 1 );
    }

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    header = (char*) malloc(50*sizeof(50));
    sprintf(header, "P5\n%d %d\n%d\n", N_x, N_y, I_max);
    header_size = strlen(header);

    job_height = world_size == 1 ? N_y : 20;
    if (world_size > 1) {
        job_remainder = (N_y % 20);
        job_remainder_size = job_remainder * N_x;
    }
    job_size = job_height * N_x;
    printf("job_remainder: %u; job_height: %u\n", job_remainder, job_height);
}

int main(int argc, char** argv) {
    init_env(argc, argv);
    init_MPI_env(argc, argv);

    pid == MASTER ? master() : slave();
    MPI_File_close(&output_file);
    printf("Process %d is terminating..\n", pid);
    MPI_Finalize();
    return 0;
}
