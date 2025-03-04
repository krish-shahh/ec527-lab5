/****************************************************************************
   Compilation Command:
   gcc -pthread -O2 -std=gnu11 test_SOR_mt.c -lm -lrt -o test_SOR_mt
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

#define CPNS 2.0    /* Cycles per nanosecond - adjust for CPU frequency */
#define GHOST 2     /* Extra rows/columns for ghost zone */
#define A 20        /* Adjusted coefficient to avoid powers of 2 */
#define B 50
#define C 70
#define NUM_TESTS 5 /* Number of different array sizes to test */
#define MAX_THREADS 8 /* Maximum number of threads */
#define MINVAL 0.0
#define MAXVAL 10.0
#define TOL 0.00001
#define OMEGA 1.75  /* Best relaxation parameter from Part 1 */

typedef double data_t;

typedef struct {
    long int rowlen;
    data_t *data;
} arr_rec, *arr_ptr;

typedef struct {
    int thread_id;
    arr_ptr v;
    int start_row;
    int end_row;
    int iterations;
} thread_data_t;

pthread_barrier_t barrier;

/* Function Prototypes */
arr_ptr new_array(long int row_len);
void init_array_rand(arr_ptr v, long int row_len);
void SOR_serial(arr_ptr v, int *iterations);
void *SOR_thread_strip(void *arg);
void *SOR_thread_interleaved(void *arg);
double interval(struct timespec start, struct timespec end);

/* Timer function */
double interval(struct timespec start, struct timespec end) {
    struct timespec temp;
    temp.tv_sec = end.tv_sec - start.tv_sec;
    temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    if (temp.tv_nsec < 0) {
        temp.tv_sec -= 1;
        temp.tv_nsec += 1000000000;
    }
    return ((double)temp.tv_sec + (double)temp.tv_nsec * 1.0e-9);
}

/* Create and initialize an array */
arr_ptr new_array(long int row_len) {
    arr_ptr result = (arr_ptr)malloc(sizeof(arr_rec));
    if (!result) return NULL;
    result->rowlen = row_len;
    result->data = (data_t *)calloc(row_len * row_len, sizeof(data_t));
    return result;
}

void init_array_rand(arr_ptr v, long int row_len) {
    srandom(row_len);
    for (long int i = 0; i < row_len * row_len; i++) {
        v->data[i] = ((double)random() / RAND_MAX) * (MAXVAL - MINVAL) + MINVAL;
    }
}

/* Standard Serial SOR */
void SOR_serial(arr_ptr v, int *iterations) {
    long int rowlen = v->rowlen;
    data_t *data = v->data;
    double change, total_change = 1.0e10;
    int iters = 0;

    while ((total_change / (rowlen * rowlen)) > TOL) {
        iters++;
        total_change = 0;
        for (long int i = 1; i < rowlen - 1; i++) {
            for (long int j = 1; j < rowlen - 1; j++) {
                change = data[i * rowlen + j] - 0.25 * (data[(i - 1) * rowlen + j] +
                                                        data[(i + 1) * rowlen + j] +
                                                        data[i * rowlen + j + 1] +
                                                        data[i * rowlen + j - 1]);
                data[i * rowlen + j] -= change * OMEGA;
                total_change += fabs(change);
            }
        }
    }
    *iterations = iters;
}

/* Strip-based Multithreaded SOR */
void *SOR_thread_strip(void *arg) {
    thread_data_t *data = (thread_data_t *)arg;
    arr_ptr v = data->v;
    long int rowlen = v->rowlen;
    double change, total_change;
    int iters = 0;

    do {
        iters++;
        total_change = 0;
        for (long int i = data->start_row; i < data->end_row; i++) {
            for (long int j = 1; j < rowlen - 1; j++) {
                change = v->data[i * rowlen + j] - 0.25 * (v->data[(i - 1) * rowlen + j] +
                                                           v->data[(i + 1) * rowlen + j] +
                                                           v->data[i * rowlen + j + 1] +
                                                           v->data[i * rowlen + j - 1]);
                v->data[i * rowlen + j] -= change * OMEGA;
                total_change += fabs(change);
            }
        }
        pthread_barrier_wait(&barrier);
    } while ((total_change / (rowlen * rowlen)) > TOL);

    data->iterations = iters;
    pthread_exit(NULL);
}

/* Interleaved Row Multithreaded SOR */
void *SOR_thread_interleaved(void *arg) {
    thread_data_t *data = (thread_data_t *)arg;
    arr_ptr v = data->v;
    long int rowlen = v->rowlen;
    double change, total_change;
    int iters = 0;

    do {
        iters++;
        total_change = 0;
        for (long int i = data->thread_id + 1; i < rowlen - 1; i += MAX_THREADS) {
            for (long int j = 1; j < rowlen - 1; j++) {
                change = v->data[i * rowlen + j] - 0.25 * (v->data[(i - 1) * rowlen + j] +
                                                           v->data[(i + 1) * rowlen + j] +
                                                           v->data[i * rowlen + j + 1] +
                                                           v->data[i * rowlen + j - 1]);
                v->data[i * rowlen + j] -= change * OMEGA;
                total_change += fabs(change);
            }
        }
        pthread_barrier_wait(&barrier);
    } while ((total_change / (rowlen * rowlen)) > TOL);

    data->iterations = iters;
    pthread_exit(NULL);
}

/* Main Function */
int main(int argc, char *argv[]) {
    struct timespec time_start, time_stop;
    double serial_time, strip_time, interleaved_time;
    int serial_iterations, strip_iterations, interleaved_iterations;

    long int array_sizes[] = {512, 2048};  // One in L3 cache, one larger than L3
    int num_threads = 4;

    for (int s = 0; s < 2; s++) {
        long int size = array_sizes[s];
        printf("\nTesting SOR on Grid Size: %ld\n", size);
        arr_ptr v0 = new_array(size);
        init_array_rand(v0, size);

        /* Serial SOR */
        clock_gettime(CLOCK_REALTIME, &time_start);
        SOR_serial(v0, &serial_iterations);
        clock_gettime(CLOCK_REALTIME, &time_stop);
        serial_time = interval(time_start, time_stop);
        printf("Serial SOR: %lf seconds, %d iterations\n", serial_time, serial_iterations);

        /* Strip-based Multithreaded SOR */
        pthread_t threads[num_threads];
        thread_data_t thread_data[num_threads];
        pthread_barrier_init(&barrier, NULL, num_threads);

        clock_gettime(CLOCK_REALTIME, &time_start);
        for (int i = 0; i < num_threads; i++) {
            thread_data[i].thread_id = i;
            thread_data[i].v = v0;
            thread_data[i].start_row = (i * size) / num_threads;
            thread_data[i].end_row = ((i + 1) * size) / num_threads;
            pthread_create(&threads[i], NULL, SOR_thread_strip, &thread_data[i]);
        }
        for (int i = 0; i < num_threads; i++) {
            pthread_join(threads[i], NULL);
        }
        clock_gettime(CLOCK_REALTIME, &time_stop);
        strip_time = interval(time_start, time_stop);
        pthread_barrier_destroy(&barrier);
        printf("Strip-Based SOR: %lf seconds\n", strip_time);

        free(v0);
    }

    return 0;
}
