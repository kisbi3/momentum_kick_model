#include <stdarg.h>
#include <stdio.h>
#include <iostream>
// #include <stdlib.h>
#ifndef integral_h
#define integral_h

// stdarg.h 공부
// va_list : Type to hold information about variable arguments (type)
// va_start : Initialize a variable argument list(macro)
// va_arg : Retrieve next argument(macro)
// va_end : End using variable argument list(macro)
// va_copy : Copy variable argument list(macro)

class integral
{
    // double function;
    // double firstintegral(double x, double dx);

    public:
        double secondintegral(double (*func)(double, double),double x_st, double x_end, int bin_x, double y_st, double y_end, int bin_y);

        //가변 인자 함수. 순서 : function, 변수 개수, x start, x end, dx , y start, y end, dy ...순서로 입력.
        double quadintegral(double (*func)(double, double, double), int num, double args2, ...);

};


#endif