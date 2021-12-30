#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <thread>
#include <time.h>

// 11/11 cuda 연습 - 여기까지.



// __global__ void function(double *x, int i, double dist){
//     dist = 3*x[i]*x[i];
// }

double function(double x, double y){
    // printf("8\n");
    // printf("%f\n",x*x+y*y);
    // *funcal = x*x+y*y;
    // double *ptr = x*x+y*y;
    return x*x+y*y;
    // return ptr;
}

//gpu를 이용하기 위해서는 '__global__'실행하려는 함수 앞에 붙여야 함.
__global__ void saxpy(int n, double dx, double dy, double *cal, double *sum){
    // int i = blockIdx.x;
    // int j = threadIdx.x;
    // printf("1");
    int index = blockIdx.x*n+threadIdx.x;
    // printf("%f\t%f\n",cal[index],sum[index]);
    sum[index] = cal[index]*dx*dy;
}

// double totalintegral = 0.;

void func(int n, double a_1, double b_1, double a_2, double b_2, double *totalintegral){
    // double a_1 = 0., b_1 = 10.;
    // double a_2 = 0., b_2 = 10.;
    // double a_3 = 0., b_3 = 5.;

    double dx, dy, x, y;

    // double *h_x = (double*)malloc(n*(sizeof(double)));
    // double *h_y = (double*)malloc(n*(sizeof(double)));
    double *h_sum = (double*)malloc(n*n*(sizeof(double)));
    double *h_cal = (double*)malloc(n*n*(sizeof(double)));

    // double x[10000] = {0.}, sum;

    // double *x;
    // double *y;
    double *sum;
    double *cal;

    // cudaMalloc((void**) &x, n*(sizeof(double)));
    // cudaMalloc((void**) &y, n*(sizeof(double)));
    cudaMalloc((void**) &sum, n*n*(sizeof(double)));
    cudaMalloc((void**) &cal, n*n*(sizeof(double)));

    // std::cout<<'1'<<std::endl;
    dx = (b_1-a_1)/double(n);
    dy = (b_2-a_2)/double(n);
    // for(int i = 0; i<n; i++){
    //     for(int j = 0; j<n; j++){
    //         h_x[i][j] = a_1+a_2 + ((b_1-a_1)/n)*i+((b_2-a_2)/n)*j;
    //     } 
    //     // std::cout<<h_x[i]<<std::endl;
    // }

    // std::cout<<dx<<std::setw(20)<<dy<<std::endl;

    // for (int i = 0; i<n; i++){
    //     h_x[i] = a_1 + ((b_1-a_1)/n)*i;
    //     h_y[i] = a_2 + ((b_2-a_2)/n)*i;

    //     // std::cout<<h_x[i]<<std::setw(20)<<h_y[i]<<std::endl;
    // }
    // printf("1\n");

    // double funcal;
    x = a_1;
    for(int i=0; i<n; i++){
        // printf("%f\n",x);
        y = a_2;
        for(int j=0; j<n; j++){
            // printf("%f\n",y);
            // funcal = function(x,y);
            // double *funcall = function(x,y,funcal);
            // printf("%f\n",funcal);
            // h_cal[i*n+j] = funcal;
            h_cal[i*n+j] = function(x,y);
            // printf("9\n");
            // printf("%f\t", h_cal[i*n+j]);
            y += dy;
        }
        // printf("\n");
        x += dx;
    }
    // printf("2\n");
    // std::cout<<std::endl<<std::endl;

    // std::cout<<'1'<<std::endl;

    // cudaMemcpy(x, h_x, n*(sizeof(double)), cudaMemcpyHostToDevice);
    // cudaMemcpy(y, h_y, n*(sizeof(double)), cudaMemcpyHostToDevice);
    cudaMemcpy(sum, h_sum, n*n*(sizeof(double)), cudaMemcpyHostToDevice);
    cudaMemcpy(cal, h_cal, n*n*(sizeof(double)), cudaMemcpyHostToDevice);

    // cal = function()

    saxpy<<<n,n>>>(n, dx, dy, cal, sum);

    cudaMemcpy(h_sum, sum, n*n*(sizeof(double)), cudaMemcpyDeviceToHost);

    double total = 0.;
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            double k = h_sum[i+n*j];
            // printf("%f",k);
            total += k;
            // std::cout<<h_sum[i]<<std::endl;            
        }
        
    }

    *totalintegral += total;

    // printf("\n%f\t%f\n",total,*totalintegral);

    // cudaFree(x);
    cudaFree(sum);
    cudaFree(cal);

    // free(h_x);
    free(h_sum);
    free(h_cal);
}

int main(void)
{
    // printf("1");
    time_t start, end;

    start = time(NULL);

    // std::cout<<'1'<<std::endl;
    int n = 1000;

    // double a_1 = 0., b_1 = 10.;
    // double a_2 = 0., b_2 = 10.;

    double *totalintegral = (double*)malloc(sizeof(double));
    // double *totalintegral;
    
    // printf("1\n");
    // func(n,0.,10.,0.,10.,totalintegral);
    // printf("1\n");
    std::thread t1(func, n, 0., 5., 0., 5., totalintegral);
    std::thread t2(func, n, 0., 5., 5., 10., totalintegral);
    std::thread t3(func, n, 5., 10., 0., 5., totalintegral);
    std::thread t4(func, n, 5., 10., 5., 10., totalintegral);

    t1.join();
    t2.join();
    t3.join();
    t4.join();

    // printf("\n1\n");

    std::cout<<*totalintegral<<std::endl;

    end = time(NULL);
    std::cout<<"걸린 시간 : "<<double(end-start)<<std::endl;

    return 0;
}