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
#include <stdbool.h>
#include <sys/syscall.h>

//mpi tags
#define _GNU_SOURCE
#define DATA 0
#define RESULT 1
#define TERMINATE 2
#define REMAINDER 3
#define MASTER 0
//uncomment this to check the load balancing of the threads
#define LOAD_BALANCE

unsigned int I_max, N_x, N_y, job_height, job_remainder, job_remainder_size, job_size, header_size;
int pid, world_size, nthreads;
double x_L, x_R, y_L, y_R;
double d_x, d_y;
char* header;
double** timer_threads;
double starting_time;
MPI_File output_file;
MPI_Status file_stat;

struct complex
{
    double real;
    double imag;
};
void init_env(int argc, char** argv);
void init_MPI_env(int argc, char** argv);
void master();
void slave();
void _worker(unsigned int start, unsigned int work_amount);
unsigned int compute_pixel(struct complex c, unsigned int max_iter); 

int main(int argc, char** argv) {
    init_env(argc, argv);
    init_MPI_env(argc, argv);
    #ifdef _OPENMP
        #pragma omp parallel
        {
            #pragma omp master
            {
                nthreads = omp_get_num_threads();
                printf("omp with %d threads\n", nthreads );
            }
            int me = omp_get_thread_num();
            #pragma omp critical
                printf("thread %2d is running on core %2d\n", me, get_cpu_id() );    
        }
        #ifdef LOAD_BALANCE
            timer_threads = (double**) malloc(sizeof(double*) * world_size);
            for (int i = 0; i < world_size; i++)
            {
                *(timer_threads + i) = (double*) malloc(sizeof(double) * nthreads);
            }
        #endif
        starting_time = omp_get_wtime();
    #else
        starting_time = MPI_Wtime();
    #endif

    pid == MASTER ? master() : slave();
    MPI_File_close(&output_file);

    #if defined(_OPENMP) && defined(LOAD_BALANCE)
        printf("Process %d; WTime: %f\n", pid, omp_get_wtime() - starting_time);
        for (int i = 0; i < nthreads; i++)
        {
            printf("PID: %d thread %d time %f\n", pid, i, timer_threads[pid][i]);
        }
        
    #else
        printf("Process %d; WTime: %f\n", pid, MPI_Wtime() - starting_time);
    #endif

    MPI_Finalize();
    return 0;
}

/**
 * Initialise enviroment
 */
void init_env(int argc, char** argv) 
{
    // resolution
    N_x = atoi(argv[1]), N_y = atoi(argv[2]);
    //area
    x_L = atof(argv[3]), y_L = atof(argv[4]);
    x_R = atof(argv[5]), y_R = atof(argv[6]);
    //maxium iterations
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
    //header for the pgm file
    header = (char*) malloc(50*sizeof(char));
    sprintf(header, "P5\n%d %d\n%d\n", N_x, N_y, I_max);
    header_size = strlen(header);

    job_height = world_size == 1 ? N_y : 20;
    if (world_size > 1) {
        job_remainder = (N_y % 20);
        job_remainder_size = job_remainder * N_x;
    }
    job_size = job_height * N_x;
}

/**
 * method that will be executed by the master only, which will act
 * as a manager that redistributes the work
 */
void master()
{
    MPI_Status stat, file_stat;
    //open the file only in the master and write the header
    MPI_File_open(MPI_COMM_SELF, "mandelbrot_set_parallel.pgm",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
    MPI_File_write(output_file, header, header_size, MPI_CHAR, &file_stat);
    MPI_File_close(&output_file);

    //collective open
    MPI_File_open(MPI_COMM_WORLD, "mandelbrot_set_parallel.pgm",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
    if (world_size == 1) {
        _worker(MASTER, N_y);
        return;
    }

    unsigned int start;
    unsigned int working_slaves = 1, row_offset = 0;
    int tag;

    //distribute the first world_size-1 jobs
    for (; working_slaves < world_size && row_offset < N_x; working_slaves++, row_offset += job_height) {
        tag = (row_offset + job_height > N_y) ? REMAINDER : DATA;
        MPI_Send(&row_offset, 1, MPI_UNSIGNED, working_slaves, tag, MPI_COMM_WORLD);
    }

    //loop in which we assign the remaining jobs to the slaves as they finish the previous assigned
    //until there are no more jobs to be assigned
    do {
        short done;
        MPI_Recv(&done, 1, MPI_SHORT, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
        int slave = stat.MPI_SOURCE;
        working_slaves--;

        //if there is still work to be made
        if (row_offset < N_y) {
            //if check == true, we are issuing the last job which will have a different
            //size because N_y is not divisible by job_height (i.e. 0 < job_remainder < 20)
            bool check = row_offset + job_height > N_y;
            tag = check ? REMAINDER : DATA;
            MPI_Send(&row_offset, 1, MPI_UNSIGNED, slave, tag, MPI_COMM_WORLD);
            row_offset += check ? job_remainder : job_height;
            working_slaves++;
        } else { 
            //otherwise we let the slave terminate
            MPI_Send(&row_offset, 1, MPI_UNSIGNED, slave, TERMINATE, MPI_COMM_WORLD);
        }
        
    } while (working_slaves > 1);
}

/**
 * method used by slave processes, that will carry out the work and 
 * reply back to the master when done
 */
void slave()
{
    MPI_Status stat;
    unsigned int row_offset;
    MPI_File_open(MPI_COMM_WORLD, "mandelbrot_set_parallel.pgm",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);

    MPI_Recv(&row_offset, 1, MPI_UNSIGNED, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);


    while (stat.MPI_TAG != TERMINATE) 
    {
        unsigned int temp_job_height = stat.MPI_TAG == DATA ? job_height : job_remainder; 
        _worker(row_offset, temp_job_height);
        short done = 1;
        MPI_Send(&done, 1, MPI_SHORT, MASTER, stat.MPI_TAG, MPI_COMM_WORLD);
        MPI_Recv(&row_offset, 1, MPI_UNSIGNED, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    }
}

/**
 * methods representing the work to be carried out by a slave
 * start represent the starting index from which the elements computed
 * will be placed in the final array of pixels
 */
void _worker(unsigned int start, unsigned int work_amount)
{
    struct complex c;
    unsigned int* buffer = (unsigned int*) malloc(work_amount * N_x * sizeof(unsigned int));
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 10) private(c) collapse(2)
    #endif
    for(unsigned int i = 0; i < work_amount; i++) {
        for(unsigned int j = 0; j < N_x; j++) {
            #if defined(_OPENMP) && defined(LOAD_BALANCE)
                double s = omp_get_wtime();
            #endif
            c.imag = y_L + (i + start) * d_y;
            c.real = x_L + j * d_x;
            *(buffer + (i * N_x + j)) = compute_pixel(c, I_max);
            #if defined(_OPENMP) && defined(LOAD_BALANCE)
                timer_threads[pid][omp_get_thread_num()] += omp_get_wtime() - s;
            #endif
        }
    }

    //write the results using offset
    MPI_File_write_at(output_file, header_size * sizeof(char) + start * N_x * sizeof(unsigned int), buffer,
                      N_y * work_amount, MPI_UNSIGNED, &file_stat);
    free(buffer);
}


/**
 * Function that given a complex number c and and integer,
 * computes if c belongs to the Mandelbrot Set and returns
 * the counter used in the loop 
 */
unsigned int compute_pixel(struct complex c, unsigned int max_iter) 
{
    unsigned int count = 0;
    struct complex z;
    z.real = 0.0; z.imag = 0.0;
    double temp;

    do 
    {
        temp = z.real * z.real - z.imag * z.imag + c.real;
        z.imag = 2 * z.real * z.imag + c.imag;
        z.real = temp;
    } while ((z.real * z.real + z.imag * z.imag < 4.0) && (count++ < max_iter));

    return count;
}

//----------------------------------------------------------------------------------------

int get_cpu_id( void )
{
#if defined(_GNU_SOURCE)                              // GNU SOURCE ------------
  
  return  sched_getcpu( );

#else

#ifdef SYS_getcpu                                     //     direct sys call ---
  
  int cpuid;
  if ( syscall( SYS_getcpu, &cpuid, NULL, NULL ) == -1 )
    return -1;
  else
    return cpuid;
  
#else      

  unsigned val;
  if ( read_proc__self_stat( CPU_ID_ENTRY_IN_PROCSTAT, &val ) == -1 )
    return -1;

  return (int)val;

#endif                                                // -----------------------
#endif

}



int read_proc__self_stat( int field, int *ret_val )
/*
  Other interesting fields:
  pid      : 0
  father   : 1
  utime    : 13
  cutime   : 14
  nthreads : 18
  rss      : 22
  cpuid    : 39
  read man /proc page for fully detailed infos
 */
{
  // not used, just mnemonic
  // char *table[ 52 ] = { [0]="pid", [1]="father", [13]="utime", [14]="cutime", [18]="nthreads", [22]="rss", [38]="cpuid"};

  *ret_val = 0;

  FILE *file = fopen( "/proc/self/stat", "r" );
  if (file == NULL )
    return -1;

  char   *line = NULL;
  int     ret;
  size_t  len;
  ret = getline( &line, &len, file );
  fclose(file);

  if( ret == -1 )
    return -1;

  char *savetoken = line;
  char *token = strtok_r( line, " ", &savetoken);
  --field;
  do { token = strtok_r( NULL, " ", &savetoken); field--; } while( field );

  *ret_val = atoi(token);

  free(line);

  return 0;
}
