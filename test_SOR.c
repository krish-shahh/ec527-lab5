/****************************************************************************
   Compilation Command:
   gcc -O1 -std=gnu11 test_SOR.c -lpthread -lrt -lm -o test_SOR
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#ifdef __APPLE__
#include "apple_pthread_barrier.h"
#endif /* __APPLE__ */

#define CPNS 2.0    /* Cycles per nanosecond - adjust for your CPU frequency */
#define GHOST 2     /* Extra rows/columns for "ghost zone" */
#define A   8       /* Coefficient of x^2 */
#define B   16      /* Coefficient of x */
#define C   32      /* Constant term */
#define NUM_TESTS 5 /* Number of different array sizes to test */
#define BLOCK_SIZE 8 /* Optimal block size determined experimentally */
#define OPTIONS 4   /* Number of SOR implementations */
#define MINVAL 0.0
#define MAXVAL 10.0
#define TOL 0.00001
#define OMEGA 1.75  /* Best performing relaxation parameter from Part 1 */

typedef double data_t;

typedef struct {
    long int rowlen;
    data_t *data;
} arr_rec, *arr_ptr;

/* Function Prototypes */
arr_ptr new_array(long int row_len);
int set_arr_rowlen(arr_ptr v, long int index);
long int get_arr_rowlen(arr_ptr v);
int init_array_rand(arr_ptr v, long int row_len);
data_t *get_array_start(arr_ptr v);
void SOR(arr_ptr v, int *iterations);
void SOR_redblack(arr_ptr v, int *iterations);
void SOR_ji(arr_ptr v, int *iterations);
void SOR_blocked(arr_ptr v, int *iterations);

double interval(struct timespec start, struct timespec end)
{
    struct timespec temp;
    temp.tv_sec = end.tv_sec - start.tv_sec;
    temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    if (temp.tv_nsec < 0) {
        temp.tv_sec -= 1;
        temp.tv_nsec += 1000000000;
    }
    return (((double)temp.tv_sec) + ((double)temp.tv_nsec) * 1.0e-9);
}

/*****************************************************************************/
int main(int argc, char *argv[])
{
    int OPTION;
    struct timespec time_start, time_stop;
    double time_stamp[OPTIONS][NUM_TESTS];
    int convergence[OPTIONS][NUM_TESTS];
    int *iterations;

    long int x, n;
    long int alloc_size = GHOST + A * (NUM_TESTS - 1) * (NUM_TESTS - 1) + B * (NUM_TESTS - 1) + C;

    printf("SOR Serial Optimizations Benchmark\n");
    printf("Using OMEGA = %0.2f\n", OMEGA);

    arr_ptr v0 = new_array(alloc_size);
    iterations = (int *)malloc(sizeof(int));

    for (OPTION = 0; OPTION < OPTIONS; OPTION++) {
        const char *option_names[] = {"Standard SOR", "Red/Black SOR", "Reversed Indices SOR", "Blocked SOR"};
        printf("\nOPTION %d: %s\n", OPTION, option_names[OPTION]);

        for (x = 0; x < NUM_TESTS && (n = A * x * x + B * x + C) <= alloc_size; x++) {
            printf("  Test %ld: Grid Size = %ld\n", x, (long)(GHOST + n));
            init_array_rand(v0, GHOST + n);
            set_arr_rowlen(v0, GHOST + n);

            clock_gettime(CLOCK_REALTIME, &time_start);
            switch (OPTION) {
                case 0: SOR(v0, iterations); break;
                case 1: SOR_redblack(v0, iterations); break;
                case 2: SOR_ji(v0, iterations); break;
                case 3: SOR_blocked(v0, iterations); break;
            }
            clock_gettime(CLOCK_REALTIME, &time_stop);

            time_stamp[OPTION][x] = interval(time_start, time_stop);
            convergence[OPTION][x] = *iterations;
        }
    }

    /* Output results */
    printf("\nFinal Results (Time in ns, Iterations to Convergence):\n");
    printf("Size, SOR Time, SOR Iters, Red/Black Time, Red/Black Iters, Reversed Time, Reversed Iters, Blocked Time, Blocked Iters\n");
    for (int i = 0; i < NUM_TESTS; i++) {
        printf("%4ld", A * i * i + B * i + C);
        for (OPTION = 0; OPTION < OPTIONS; OPTION++) {
            printf(", %10.4g", (double)CPNS * 1.0e9 * time_stamp[OPTION][i]);
            printf(", %4d", convergence[OPTION][i]);
        }
        printf("\n");
    }

    free(iterations);
    return 0;
}

/*********************************/

/* Function Definitions */
arr_ptr new_array(long int row_len)
{
    arr_ptr result = (arr_ptr)malloc(sizeof(arr_rec));
    if (!result) return NULL;
    result->rowlen = row_len;
    result->data = (data_t *)calloc(row_len * row_len, sizeof(data_t));
    if (!result->data) {
        free(result);
        return NULL;
    }
    return result;
}

int set_arr_rowlen(arr_ptr v, long int row_len) { v->rowlen = row_len; return 1; }
long int get_arr_rowlen(arr_ptr v) { return v->rowlen; }
data_t *get_array_start(arr_ptr v) { return v->data; }

int init_array_rand(arr_ptr v, long int row_len)
{
    srandom(row_len);
    for (long int i = 0; i < row_len * row_len; i++) {
        v->data[i] = ((double)random() / RAND_MAX) * (MAXVAL - MINVAL) + MINVAL;
    }
    return 1;
}

/************************************/

/* Standard SOR */
void SOR(arr_ptr v, int *iterations) {
    long int rowlen = get_arr_rowlen(v);
    data_t *data = get_array_start(v);
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

/* SOR red/black */
void SOR_redblack(arr_ptr v, int *iterations)
{
  int i, j, redblack;
  long int ti;
  long int rowlen = get_arr_rowlen(v);
  data_t *data = get_array_start(v);
  double change, total_change = 1.0e10;   /* start w/ something big */
  int iters = 0;

  ti = 0;
  redblack = 0;
  /* The while condition here tests the tolerance limit *only* when
     redblack is 0, which ensures we exit only after having done a
     full update (red + black) */
  while ((redblack == 1)
        || ((total_change/(double)(rowlen*rowlen)) > (double)TOL) )
  {
    /* Reset sum of total change only when starting a black scan. */
    if (redblack == 0) {
      total_change = 0;
    }
    for (i = 1; i < rowlen-1; i++) {
      /* The j loop needs to start at j=1 on row 0 and all even rows,
         and start at j=2 on odd rows; but when redblack is true it does
         just the opposite; and it always increments by 2. */
      for (j = 1 + ((i^redblack)&1); j < rowlen-1; j+=2) {
        change = data[i*rowlen+j] - .25 * (data[(i-1)*rowlen+j] +
                                          data[(i+1)*rowlen+j] +
                                          data[i*rowlen+j+1] +
                                          data[i*rowlen+j-1]);
        data[i*rowlen+j] -= change * OMEGA;
        if (change < 0) {
          change = -change;
        }
        total_change += change;
        ti++;
      }
    }
    if (abs(data[(rowlen-2)*(rowlen-2)]) > 10.0*(MAXVAL - MINVAL)) {
      printf("SOR: SUSPECT DIVERGENCE iter = %ld\n", iters);
      break;
    }
    redblack ^= 1;
    iters++;
  }
  /* A "red scan" only updates half of the array, and likewise for a
     "black scan"; so we need to divide iters by 2 to convert our count of
     "reds+blacks" to a count of "full scans" */
  iters /= 2;
  *iterations = iters;
  printf("    SOR_redblack() done after %d iters\n", iters);
  /* printf("ti == %ld, per iter %ld\n", ti, ti/iters); */
} /* End of SOR_redblack */

/* SOR with reversed indices */
void SOR_ji(arr_ptr v, int *iterations)
{
  long int i, j;
  long int rowlen = get_arr_rowlen(v);
  data_t *data = get_array_start(v);
  double change, total_change = 1.0e10;   /* start w/ something big */
  int iters = 0;

  while ((total_change/(double)(rowlen*rowlen)) > (double)TOL) {
    iters++;
    total_change = 0;
    for (j = 1; j < rowlen-1; j++) {
      for (i = 1; i < rowlen-1; i++) {
        change = data[i*rowlen+j] - .25 * (data[(i-1)*rowlen+j] +
                                          data[(i+1)*rowlen+j] +
                                          data[i*rowlen+j+1] +
                                          data[i*rowlen+j-1]);
        data[i*rowlen+j] -= change * OMEGA;
        if (change < 0){
          change = -change;
        }
        total_change += change;
      }
    }
    if (abs(data[(rowlen-2)*(rowlen-2)]) > 10.0*(MAXVAL - MINVAL)) {
      printf("SOR_ji: SUSPECT DIVERGENCE iter = %d\n", iters);
      break;
    }
  }
  *iterations = iters;
  printf("    SOR_ji() done after %d iters\n", iters);
}

/* SOR w/ blocking */
void SOR_blocked(arr_ptr v, int *iterations)
{
  long int i, j, ii, jj;
  long int rowlen = get_arr_rowlen(v);
  data_t *data = get_array_start(v);
  double change, total_change = 1.0e10;
  int iters = 0;
  int k;

  if ((rowlen-2) % (BLOCK_SIZE)) {
    fprintf(stderr,
"SOR_blocked: Total array size must be 2 more than a multiple of BLOCK_SIZE\n"
"(because the top/right/left/bottom rows are not scanned)\n"
"Make sure all coefficients A, B, and C are multiples of %d\n", BLOCK_SIZE);
    exit(-1);
  }

  while ((total_change/(double)(rowlen*rowlen)) > (double)TOL) {
    iters++;
    total_change = 0;
    for (ii = 1; ii < rowlen-1; ii+=BLOCK_SIZE) {
      for (jj = 1; jj < rowlen-1; jj+=BLOCK_SIZE) {
        for (i = ii; i < ii+BLOCK_SIZE; i++) {
          for (j = jj; j < jj+BLOCK_SIZE; j++) {
            change = data[i*rowlen+j] - .25 * (data[(i-1)*rowlen+j] +
                                              data[(i+1)*rowlen+j] +
                                              data[i*rowlen+j+1] +
                                              data[i*rowlen+j-1]);
            data[i*rowlen+j] -= change * OMEGA;
            if (change < 0){
              change = -change;
            }
            total_change += change;
          }
        }
      }
    }
    if (abs(data[(rowlen-2)*(rowlen-2)]) > 10.0*(MAXVAL - MINVAL)) {
      printf("SOR_blocked: SUSPECT DIVERGENCE iter = %d\n", iters);
      break;
    }
  }
  *iterations = iters;
  printf("    SOR_blocked() done after %d iters\n", iters);
} /* End of SOR_blocked */

