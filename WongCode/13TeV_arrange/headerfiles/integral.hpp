#ifndef integral_h
#define integral_h


class integral
{
    // double function;
    // double firstintegral(double x, double dx);

    public:
        double secondintegral(double (*func)(double, double),double x_st, double x_end, int bin_x, double y_st, double y_end, int bin_y);
};

// class cuda_integral
// {
//     public:
//         double secondintegral(double (*func)(double, double),double x_st, double x_end, int bin_x, double y_st, double y_end, int bin_y);

//     private:
//         // int arraySize;
//         // int blockSize;
//         // void integrate(double x_st, double x_end, double *x, double y_st, double y_end, double *y, double *sum);
//         // void sumCommMultiBlock(double *gArr, int arraySize, double *gOut);  //앞에 __global__을 넣어야 하려나?
//         double twodimsumarray(double *gArr, int bin_x, int bin_y, double *gOut);
// };


#endif