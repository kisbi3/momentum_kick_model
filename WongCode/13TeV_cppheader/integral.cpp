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

// void *safe_malloc(size_t n)
// {
//     void *p = malloc(n);
//     if (p == NULL) {
//         fprintf(stderr, "Fatal: failed to allocate %zu bytes.\n", n);
//         abort();
//     }
//     return p;
// }

double integral::quadintegral(double (*func)(double , double, double), int num, double args2,...){
    // printf("\n");

    va_list ap;
    // va_start(ap, args2);    //args2의 첫 숫자는 여기에서 사라진다. args2로 입력됨.
    // double *xyz = (double*)malloc(num*sizeof(double));
    double *xyz;
    int *bin;
    int loop = 1; // loop를 돌 횟수(적분 횟수)
    bin = (int *)malloc(num*sizeof(int));
    xyz = (double *)malloc(num*sizeof(double));

    // xyz[0]=xyz[1]=xyz[2]=xyz[3]=xyz[4]=xyz[5]=xyz[6]=xyz[7]=xyz[8] = 0.;
    // bin[0]=bin[1]=bin[2]=0;
    // printf("%f\n", args2);
    for(int i=0;i<num;i++){   //args2의 첫 인자의 개수만큼 반복
        // printf("1\n");
        double x_st = 0.;
        if(i==0){
            va_start(ap, args2);
            x_st = args2;
        }
        else{
            x_st = va_arg(ap, double);
        }
        // double x_st = va_arg(ap, double);   //va_arg는 double만큼 할당된 인수를 내뱉고 바로 다음으로 넘어감.
        double x_end = va_arg(ap, double);
        double dx = va_arg(ap, double);
        bin[i] = int ((x_end-x_st)/dx);
        xyz[3*i] = x_st;
        xyz[3*i+1] = x_end;
        xyz[3*i+2] = dx;


        // printf("%d\n", i);
        // printf("%f\t%f\t%f\t%d\n", xyz[3*i], xyz[3*i+1], xyz[3*i+2], bin[i]);
        // printf("%f\t%f\t%f\n", xyz[3*i], xyz[3*i+1], xyz[3*i+2]);
        // printf("%f\t%f\t%f\t%d\n", x_st, x_end, dx, bin[i]);
        loop *= bin[i];
    }

    // printf("%d", num);
    double sum = 0.;
    //어떻게 하지?
    if(num==1){
        for(int j = 0;j<bin[0];j++){
            sum += func(xyz[0]+xyz[2]*j)*xyz[2];
        }
        free(xyz);
        free(bin);
        return sum;        
    }
    else if(num==2){
        for(int j = 0;j<bin[0];j++){
            for(int k = 0; k<bin[1];k++){
                sum += func(xyz[0]+xyz[2]*j, xyz[3]+xyz[5]*k)*xyz[2]*xyz[5];
            }
        }
        free(xyz);
        free(bin);
        return sum;        
    }


    else if(num==3){
        for(int j = 0;j<bin[0];j++){
            for(int k = 0; k<bin[1];k++){
                for(int l = 0; l<bin[2];l++){
                    sum += func(xyz[0]+xyz[2]*j, xyz[3]+xyz[5]*k, xyz[6]+xyz[8]*l)*xyz[2]*xyz[5]*xyz[8];
                }
            }
        }
        free(xyz);
        free(bin);
        return sum;        
    }


    else{
        printf("Error!\n");
        return 1;
    }

    // return xyz[4];


}