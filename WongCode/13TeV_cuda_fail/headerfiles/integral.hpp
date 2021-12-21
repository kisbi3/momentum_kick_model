#include <stdarg.h>
#include <stdio.h>
#ifndef integral_h
#define integral_h


class integral
{
    // double function;
    // double firstintegral(double x, double dx);

    public:
        double secondintegral(double (*func)(double, double),double x_st, double x_end, int bin_x, double y_st, double y_end, int bin_y);

        //가변 인자 함수. 순서 : function, x start, x end, dx , y start, y end, dy ...순서로 입력.
        double quadintegral(double (*func)(double, ...), double args, ...);
};


#endif