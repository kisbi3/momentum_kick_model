#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <thread>
#include <time.h>

// 11/11 cuda 연습 - 여기까지.



// __global__ void function(double *x, int i, double dist){
//     dist = 3*x[i]*x[i];
// }

// double function(){

// }

//gpu를 이용하기 위해서는 '__global__'실행하려는 함수 앞에 붙여야 함.
__global__ void saxpy(int n, double *x, double *y, double dx, double dy, double *sum){
    // int i = blockIdx.x*blockDim.x;

    // printf("%d", blockIdx.x);


    int i = blockIdx.x;
    int j = threadIdx.x;
    // double *func;


    // printf("%d\n",j);

    int index = blockIdx.x*n+threadIdx.x;


    // __shared__ void function(double *x, double *y, double *func){
    //     *func = x*x+y*y;
    // }
    // double dx = (b-a)/n;
    sum[index] = (x[i]*x[i]+y[j]*y[j])*dx*dy;

    // func = function(&x[i], &y[j], &func);
    // sum[index] = dist*dx*dy;

    // double x_temp = x[i];
    // double y_temp = y[j];
    // sum[index] = x_temp*y_temp;
    // sum[index] = x[i]*dx+y[j]*dy;


    // printf("%f\t%f\t%f\t%d\n",x[i],y[j],sum[index],index);


    // sum[index] = x[i]*dx;
    
    // double dist = function(x, i, dist);
    // sum[i] = dist*dx;

    // if( i < n )
    //     sum[i] = a*x[i]+b;
}

double totalintegral = 0.;

void func(int n, double a_1, double b_1, double a_2, double b_2){
    // double a_1 = 0., b_1 = 10.;
    // double a_2 = 0., b_2 = 10.;
    // double a_3 = 0., b_3 = 5.;

    double dx, dy;

    double *h_x = (double*)malloc(n*(sizeof(double)));
    double *h_y = (double*)malloc(n*(sizeof(double)));
    double *h_sum = (double*)malloc(n*n*(sizeof(double)));

    // double x[10000] = {0.}, sum;

    double *x;
    double *y;
    double *sum;

    cudaMalloc((void**) &x, n*(sizeof(double)));
    cudaMalloc((void**) &y, n*(sizeof(double)));
    cudaMalloc((void**) &sum, n*n*(sizeof(double)));

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

    for (int i = 0; i<n; i++){
        h_x[i] = a_1 + ((b_1-a_1)/n)*i;
        h_y[i] = a_2 + ((b_2-a_2)/n)*i;

        // std::cout<<h_x[i]<<std::setw(20)<<h_y[i]<<std::endl;
    }
    
    // std::cout<<std::endl<<std::endl;

    // std::cout<<'1'<<std::endl;

    cudaMemcpy(x, h_x, n*(sizeof(double)), cudaMemcpyHostToDevice);
    cudaMemcpy(y, h_y, n*(sizeof(double)), cudaMemcpyHostToDevice);
    cudaMemcpy(sum, h_sum, n*n*(sizeof(double)), cudaMemcpyHostToDevice);

    saxpy<<<n,n>>>(n, x, y, dx, dy, sum);

    cudaMemcpy(h_sum, sum, n*n*(sizeof(double)), cudaMemcpyDeviceToHost);

    double total = 0.;
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            double k = h_sum[i+n*j];
            total += k;
            // std::cout<<h_sum[i]<<std::endl;            
        }
        
    }

    totalintegral += total;

    cudaFree(x);
    cudaFree(sum);

    free(h_x);
    free(h_sum);
}

int main(void)
{
    // using std::thread;
    // using std::cout;
    // using std::endl;
    // using std::setw;
    time_t start, end;

    start = time(NULL);

    // std::cout<<'1'<<std::endl;
    int n = 10000;

    // double a_1 = 0., b_1 = 10.;
    // double a_2 = 0., b_2 = 10.;

    std::thread t1(func, n, 0., 5., 0., 5.);
    std::thread t2(func, n, 0., 5., 5., 10.);
    std::thread t3(func, n, 5., 10., 0., 5.);
    std::thread t4(func, n, 5., 10., 5., 10.);

    t1.join();
    t2.join();
    t3.join();
    t4.join();
    // int size = n*n*(sizeof(double));

    // double dx, dy;
    // // double *x, *sum;

    // double *h_x = (double*)malloc(n*(sizeof(double)));
    // double *h_y = (double*)malloc(n*(sizeof(double)));
    // double *h_sum = (double*)malloc(n*n*(sizeof(double)));

    // // double x[10000] = {0.}, sum;

    // double *x;
    // double *y;
    // double *sum;

    // cudaMalloc((void**) &x, n*(sizeof(double)));
    // cudaMalloc((void**) &y, n*(sizeof(double)));
    // cudaMalloc((void**) &sum, n*n*(sizeof(double)));

    // double a_1 = 0., b_1 = 10.;
    // double a_2 = 0., b_2 = 10.;
    // // double a_3 = 0., b_3 = 5.;

    // // std::cout<<'1'<<std::endl;
    // dx = (b_1-a_1)/double(n);
    // dy = (b_2-a_2)/double(n);
    // // for(int i = 0; i<n; i++){
    // //     for(int j = 0; j<n; j++){
    // //         h_x[i][j] = a_1+a_2 + ((b_1-a_1)/n)*i+((b_2-a_2)/n)*j;
    // //     } 
    // //     // std::cout<<h_x[i]<<std::endl;
    // // }

    // // std::cout<<dx<<std::setw(20)<<dy<<std::endl;

    // for (int i = 0; i<n; i++){
    //     h_x[i] = a_1 + ((b_1-a_1)/n)*i;
    //     h_y[i] = a_2 + ((b_2-a_2)/n)*i;

    //     // std::cout<<h_x[i]<<std::setw(20)<<h_y[i]<<std::endl;
    // }
    
    // // std::cout<<std::endl<<std::endl;

    // // std::cout<<'1'<<std::endl;

    // cudaMemcpy(x, h_x, n*(sizeof(double)), cudaMemcpyHostToDevice);
    // cudaMemcpy(y, h_y, n*(sizeof(double)), cudaMemcpyHostToDevice);
    // cudaMemcpy(sum, h_sum, n*n*(sizeof(double)), cudaMemcpyHostToDevice);

    // // std::cout<<'1'<<std::endl;



    // // saxpy<<<n,n>>>(n, x, y, dx, dy, sum);

    // // std::cout<<'1'<<std::endl;
    // // 'saxpy'의 함수를 nx3개의 gpu thread가 실행한다는 의미
    // // saxpy<<<b, n>>>();
    // // 에서 saxpy는 함수의 이름, b는 함수를 수행할 block의 개수, n은 다시 하나의 thread block 안에 몇 개의 thread가 존재하는지를 정하는 것.

    // cudaMemcpy(h_sum, sum, n*n*(sizeof(double)), cudaMemcpyDeviceToHost);

    // total = 0.;
    // for(int i = 0; i<n; i++){
    //     for(int j = 0; j<n; j++){
    //         double k = h_sum[i+n*j];
    //         total += k;
    //         // std::cout<<h_sum[i]<<std::endl;            
    //     }
        
    // }

    // cudaFree(x);
    // cudaFree(sum);

    // free(h_x);
    // free(h_sum);

    // std::cout<<total<<std::endl;
    std::cout<<totalintegral<<std::endl;

    end = time(NULL);
    // std::cout<<double(end-start)<<std::endl;
    std::cout<<"걸린 시간 : "<<double(end-start)<<std::endl;

    return 0;
}