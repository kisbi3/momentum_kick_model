#include <stdio.h>
#ifndef integral_cuh
#define integral_cuh
// #include <cuda_runtime.h>

// extern "C" double cuda_secondintegral(double (*func)(double, double), double x_st, double x_end, int bin_x, double y_st, double y_end, int bin_y);
extern "C" double cuda_secondintegral(double x_st, double x_end, int bin_x, double y_st, double y_end, int bin_y);
extern "C" double cuda_twodimsumarray(double *gArr, int bin_x, int bin_y, double *gOut);
// extern "C" double cuda_function(double (*func)(double, double), double x, double y);

// extern "C" double cuda_func(double x, double y);



#endif