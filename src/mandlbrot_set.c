#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
//#include <complex.h>

//mpi tags
#define DATA 0
#define RESULT 1
#define TERMINATE 2
#define MASTER 0
//constant for 
#define LOAD_BALANCE_TIME true

int N_x, N_y, I_max, job_width, job_size, *image, *result;
int pid, world_size;
double x_L, x_R, y_L, y_R;
double d_x, d_y;

struct complex
{
    double real;
    double imag;
};

int cal_pixel(struct complex c, int max_iter) {
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



void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name)
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

void _worker(int start, int* result)
{
    //int *color = result + 1;
    struct complex c;
    result[0] = start;
    #pragma omp parallel for schedule(dynamic, 10) private(c) collapse(2)
    for(int j=0; j < N_y; j++) {
        for(int i = 0; i < job_width; i++) {
            #if LOAD_BALANCE_TIME
                double t = omp_get_wtime();
            #endif
            c.real = i * d_x + x_L;
            c.imag = j * d_y + y_L;
            result[j * job_width + i + 1] = cal_pixel(c, I_max);
            #if LOAD_BALANCE_TIME
                thread_time[pid][omp_get_thread_num()] = omp_get_wtime() - t;
            #endif
        }
    }
    
    /*for (int i = start; i < job_width*start; i++) {
        for (int y = 0; y < height; y++) {
            #ifdef _LODE_BALANCE_ANALYSIS_
                double s = omp_get_wtime();
            #endif
                c.real = i * d_x + x_L;
                c.imag = y * d_y + y_L;
                result[y * job_width + i] = calc_pixel(c);
            #ifdef _LODE_BALANCE_ANALYSIS_
                timer[pid][omp_get_thread_num()] += omp_get_wtime() - s;
            #endif
        }
    }*/
    //result[0] = start;
}

void master()
{
    //if (gui) create_display(0, 0, height, width);
    if (world_size == 1) {
        _worker(MASTER, result);
        //if (gui) { gui_draw(result[0], result + 1); flush(); }
        return;
    }

    MPI_Status stat;
    int actives = 1, jobs = 0, start;
    for (; actives < world_size && jobs < N_x; actives++, jobs += job_width) {
        MPI_Send(&jobs, 1, MPI_INT, actives, DATA, MPI_COMM_WORLD);
    }
    
    do {
        MPI_Recv(result, job_size+1, MPI_INT, MPI_ANY_SOURCE, RESULT, MPI_COMM_WORLD, &stat);
        int slave = stat.MPI_SOURCE;
        start = result[0];
        //no need to parallelize this loop, just 20 iterations
        for(int i = 1; start < start + job_width; i++, start++) {
            image[start] = result[i];
        }
        
        actives--;
        if (jobs < N_x) {
            MPI_Send(&jobs, 1, MPI_INT, slave, DATA, MPI_COMM_WORLD);
            jobs += job_width;
            actives++;
        } else { 
            MPI_Send(&jobs, 1, MPI_INT, slave, TERMINATE, MPI_COMM_WORLD);
        }
    } while (actives > 1);
    write_pgm_image((void*) image, I_max, N_x, N_y, "mandelbrot_set")
}

void slave()
{
    MPI_Status stat;
    unsigned int col;
    MPI_Recv(col, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    while (stat.MPI_TAG == DATA) {
        _worker(col, result);
        MPI_Send(result, job_size+1, MPI_INT, MASTER, RESULT, MPI_COMM_WORLD);
        MPI_Recv(col, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    }
}

/**
 * Initialise enviroment
 **/
void initial_env(int argc, char** argv) {
    N_x = stoi(argv[1]), N_y = stoi(argv[2]);
    x_L = stod(argv[3]), y_L = stod(argv[4]);
    x_R = stod(argv[5]), y_R = stod(argv[6]);
    I_max = stoi(argv[7]);

    d_x = (x_R - x_L) / N_x;
    d_y = (y_R - y_L) / N_y;
}

/**
 * Initialise MPI environment
 **/
void initial_MPI_env(int argc, char** argv)
{
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    if(pid == MASTER) {
        image = (int*) malloc(N_y * N_x * sizeof(int));
    }
    
    job_width = world_size == 1 ? N_x : 20;
    job_size = job_width * N_y;
    result = (int*) malloc(job_size + 1 * sizeof(int));

    #if LOAD_BALANCE_TIME
        //this matrix is 20x20, no need to allocate it in the heap
        int thread_time[world_size][omp_get_num_threads];
    #endif
}

void start()
{
    pid == MASTER ? master() : slave();
}

int main(int argc, char** argv) {
    initial_env(argc, argv);
    initial_MPI_env(argc, argv);
    start();
    #if LOAD_BALANCE_TIME
        for(int i = 0; i < 20; i++) {
            printf("Thread times for process %d are: \n", i);
            for(int j = 0; j < 20; j++) {
                printf("thread num %d: %d", j, thread_time[i][j]);
            }
        }
    #endif
}