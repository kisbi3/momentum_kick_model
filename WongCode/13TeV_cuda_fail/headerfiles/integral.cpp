#include "integral.hpp"

// double integral::firstintegral(double x, double dx){

// }

double integral::secondintegral(double (*func)(double, double), double x_st, double x_end, int bin_x, double y_st, double y_end, int bin_y){
    int i, j;
    double x, y, sum = 0.;
    double dx = (x_end-x_st)/bin_x;
    double dy = (y_end-y_st)/bin_y;
    for(i=0;i<bin_x;i++){
        x = x_st+dx*i;
        for(int j=0;j<bin_y;j++){
            y = y_st+dy*j;
            sum += func(x, y)*dx*dy;
        }
        y = y_st;
    }

    return sum;
}

double integral::quadintegral(double (*func)(double, ...), double args, ...){
    va_list ap;

    va_start(ap, args);
    for(int i=0;i<args;i++){
        int num = va_arg(ap, int);  
        printf("%d", num);
    }


    va_end(ap);
}