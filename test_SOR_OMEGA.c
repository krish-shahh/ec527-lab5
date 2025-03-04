/*****************************************************************************

   gcc -O1 test_SOR_OMEGA.c -lm -o test_SOR_OMEGA

 */

 #include <limits.h>
 #include <math.h>
 #include <pthread.h>
 #include <stdio.h>
 #include <stdlib.h>
 
 #define MINVAL   0.0
 #define MAXVAL  100.0
 
 #define TOL 0.00001
 
 #define START_OMEGA 0.50 /* The first OMEGA value to try */
 #define OMEGA_INC 0.01   /* OMEGA increment for each O_ITERS */
 #define O_ITERS 150      /* How many OMEGA values to test */
 
 #define PER_O_TRIALS 10  /* trials per OMEGA value */
 
 typedef double data_t;
 
 /* Create abstract data type for a 2D array */
 typedef struct {
     long int rowlen;
     data_t *data;
 } arr_rec, *arr_ptr;
 
 /* Function Prototypes */
 arr_ptr new_array(long int row_len);
 int set_arr_rowlen(arr_ptr v, long int index);
 long int get_arr_rowlen(arr_ptr v);
 int init_array(arr_ptr v, long int row_len);
 int init_array_rand(arr_ptr v, long int row_len);
 int print_array(arr_ptr v);
 double fRand(double fMin, double fMax);
 void SOR(arr_ptr v, int *iterations);
 
 /* Define different array sizes for testing */
 #define NUM_ARRAY_SIZES 4
 int array_sizes[NUM_ARRAY_SIZES] = {32, 64, 128, 256}; // Small to large arrays
 
 double OMEGA;  // Declare OMEGA globally
 
 /*****************************************************************************/
 int main(int argc, char *argv[])
 {
     double convergence[O_ITERS][NUM_ARRAY_SIZES];  
     int *iterations;
     long int i, j, k;
 
     printf("SOR OMEGA test\n");
 
     /* Allocate memory for iterations count */
     iterations = (int *) malloc(sizeof(int));
 
     for (k = 0; k < NUM_ARRAY_SIZES; k++) {
         int current_size = array_sizes[k];
         printf("\nTesting Array Size: %dx%d\n", current_size, current_size);
         
         /* Declare and initialize the array */
         arr_ptr v0 = new_array(current_size);
         if (!v0) {
             printf("Memory allocation failed for array size %d\n", current_size);
             exit(EXIT_FAILURE);
         }
 
         OMEGA = START_OMEGA;
         for (i = 0; i < O_ITERS; i++) {
             printf("%0.2f", OMEGA);
             double acc = 0.0;
             for (j = 0; j < PER_O_TRIALS; j++) {
                 init_array_rand(v0, current_size);
                 SOR(v0, iterations);
                 acc += (double)(*iterations);
                 printf(", %d", *iterations);
             }
             printf("\n");
             convergence[i][k] = acc / (double)(PER_O_TRIALS);
             OMEGA += OMEGA_INC;
         }
 
         free(v0->data);
         free(v0);
     }
 
     /* Print results for graphing */
     printf("\nOMEGA, ");
     for (k = 0; k < NUM_ARRAY_SIZES; k++)
         printf("Array Size %dx%d, ", array_sizes[k], array_sizes[k]);
     printf("\n");
 
     OMEGA = START_OMEGA;
     for (i = 0; i < O_ITERS; i++) {
         printf("%0.4f", OMEGA);
         for (k = 0; k < NUM_ARRAY_SIZES; k++)
             printf(", %0.1f", convergence[i][k]);
         printf("\n");
         OMEGA += OMEGA_INC;
     }
 
     free(iterations);
     return 0;
 }
 
 /*********************************/
 
 /* Create 2D array of specified length per dimension */
 arr_ptr new_array(long int row_len)
 {
     arr_ptr result = (arr_ptr) malloc(sizeof(arr_rec));
     if (!result) {
         return NULL;
     }
     result->rowlen = row_len;
 
     if (row_len > 0) {
         data_t *data = (data_t *) calloc(row_len * row_len, sizeof(data_t));
         if (!data) {
             free(result);
             return NULL;
         }
         result->data = data;
     } else {
         result->data = NULL;
     }
     return result;
 }
 
 /* Set row length of array */
 int set_arr_rowlen(arr_ptr v, long int row_len)
 {
     v->rowlen = row_len;
     return 1;
 }
 
 /* Return row length of array */
 long int get_arr_rowlen(arr_ptr v)
 {
     return v->rowlen;
 }
 
 /* initialize 2D array with incrementing values (0.0, 1.0, 2.0, 3.0, ...) */
 int init_array(arr_ptr v, long int row_len)
 {
     long int i;
 
     if (row_len > 0) {
         v->rowlen = row_len;
         for (i = 0; i < row_len * row_len; i++) {
             v->data[i] = (data_t)(i);
         }
         return 1;
     }
     return 0;
 }
 
 /* initialize array with random numbers */
 int init_array_rand(arr_ptr v, long int row_len)
 {
     long int i;
     if (row_len > 0) {
         v->rowlen = row_len;
         for (i = 0; i < row_len * row_len; i++) {
             v->data[i] = (data_t)(fRand((double)(MINVAL), (double)(MAXVAL)));
         }
         return 1;
     }
     return 0;
 }
 
 /* Generate a random double in range */
 double fRand(double fMin, double fMax)
 {
     double f = (double)rand() / RAND_MAX;
     return fMin + f * (fMax - fMin);
 }
 
 /* print all elements of an array */
 int print_array(arr_ptr v)
 {
     long int i, j, row_len = v->rowlen;
 
     printf("row length = %ld\n", row_len);
     for (i = 0; i < row_len; i++) {
         for (j = 0; j < row_len; j++)
             printf("%.4f ", (data_t)(v->data[i * row_len + j]));
         printf("\n");
     }
     return 0;
 }
 
 /************************************/
 
 /* SOR */
 void SOR(arr_ptr v, int *iterations)
 {
     long int i, j;
     long int row_len = get_arr_rowlen(v);
     data_t *data = v->data;
     double change, mean_change = 1.0e10;
     int iters = 0;
 
     if (OMEGA >= 2.0 || OMEGA < 0.1) {
         *iterations = INT_MAX;
         return;
     }
 
     while ((mean_change / (double)(row_len * row_len)) > (double)TOL) {
         iters++;
         mean_change = 0;
         for (i = 1; i < row_len - 1; i++) {
             for (j = 1; j < row_len - 1; j++) {
                 change = data[i * row_len + j]
                     - 0.25 * (data[(i - 1) * row_len + j] +
                               data[(i + 1) * row_len + j] +
                               data[i * row_len + j + 1] +
                               data[i * row_len + j - 1]);
                 data[i * row_len + j] -= change * OMEGA;
                 mean_change += fabs(change);
             }
         }
     }
     *iterations = iters;
 }
 