#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "openacc.h"
#define SIZE 3000
//#define VEC_L 2000

void Matrix_Mul_p ( float *restrict c, float* a, float* b, int N ) {
#pragma acc kernels present (c, a, b)
    for ( int n = 0; n < N; n++ ) {
       #pragma acc loop
        for ( int m = 0; m < N; m++ ) {
         #pragma acc loop   //float sum = 0.f;
         for ( int k = 0; k < N; k++)  c[m + n * N] += a[k + n * N] * b[k + m * N];
        }
    }
}

void Matrix_Mul ( float *restrict c, float* a, float* b, int N ) {

    for ( int n = 0; n < N; n++ ) {
       
        for ( int m = 0; m < N; m++ ) {
           //float sum = 0.f;
         for ( int k = 0; k < N; k++)  c[m + n * N] += a[k + n * N] * b[k + m * N];
        }
    }
}



int main ()
{
    float* a, *b, *c;
    clock_t t;
    a = (float*) malloc(SIZE*SIZE*sizeof(float));
    b = (float*) malloc(SIZE*SIZE*sizeof(float));
    c = (float*) malloc(SIZE*SIZE*sizeof(float));
    for ( int i = 0; i < SIZE*SIZE; i++ ) a[i] = 1.;
    for ( int j = 0; j < SIZE*SIZE; j++ ) b[j] = 0.1;
    for ( int j = 0; j < SIZE*SIZE; j++ ) c[j] = 0;
     if(SIZE <= 3000) {
     t = clock();
    Matrix_Mul (c, a, b, SIZE);
    t = clock()-t;
    printf("%lf, without parallelazation time: %lf sec\n",c[0],t/(float)CLOCKS_PER_SEC);
     }
     for ( int j = 0; j < SIZE*SIZE; j++ ) c[j] = 0;
    t = clock();
 #pragma acc data copyin (a[0:SIZE*SIZE], b[0:SIZE*SIZE]) copyout (c[0:SIZE*SIZE])
 {
   Matrix_Mul_p (c, a, b, SIZE);
 }
    t = clock()-t;
    printf("%lf, time: %lf sec\n", c[0], t/(float)CLOCKS_PER_SEC);
    //printf("%s",acc_get_device_type());
 return 0;
}




