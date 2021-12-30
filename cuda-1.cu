#include <iostream>

// __global__ void function(float *x, int i, float dist){
//     dist = 3*x[i]*x[i];
// }

__device__ double function2(float x){
    double sqrSnn = 200.;
    double mp = 0.938272046;
    double m = 0.13957018;

    // double x = rapiditxnit(x);
    double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mpi
    // printf("2");
    double squareroot=sqrt(m*m+x*x);
    // printf("2");
    double xabs = fabs(x);
    double result = (squareroot/m)*exp(xabs-yb);

    return result;
}

__device__ double functions(float x){
    double y = function2(x)*x;
    double a = 0.5;
    double T = 0.5;
    double md = 1.;
    double m = 0.13957018;
    // return y;
    if(y>=1.){
        return 0.;
        // cal[index] = 0.;
    }
    else{
        return x*pow(1-y,a)*exp(-sqrt(m*m+x*x)/T)/sqrt(md*md+x*x);
        // cal[index] = x*pow(1-x,a)*exp(-sqrt(m*m+x*x)/T)/sqrt(md*md+x*x);
    }
}


//gpu를 이용하기 위해서는 '__global__'실행하려는 함수 앞에 붙여야 함.
__global__ void saxpy(int n, float a, float b, float *x, float dx, float *sum){
    // int i = blockIdx.x*blockDim.x;
    int i = blockIdx.x;
    // float dx = (b-a)/n;

    double fx = functions(x[i]);
    // printf("%f\n", fx);
    sum[i] = (fx)*dx;
    // sum[i] = exp(x[i])*dx;
    
    // float dist = function(x, i, dist);
    // sum[i] = dist*dx;

    // if( i < n )
    //     sum[i] = a*x[i]+b;
}

int main(void)
{
    // std::cout<<'1'<<std::endl;
    int n = 1000000;
    float dx;
    float total;
    // float *x, *sum;

    float *h_x = (float*)malloc(n*sizeof(float));
    float *h_sum = (float*)malloc(n*sizeof(float));

    // float x[10000] = {0.}, sum;

    float *x;
    float *sum;

    cudaMalloc((void**) &x, n*sizeof(float));
    cudaMalloc((void**) &sum, n*sizeof(float));

    float a = 0., b = 5.;

    // std::cout<<'1'<<std::endl;
    dx = (b-a)/float(n);
    for(int i = 0; i<n; i++){
        h_x[i] = a + ((b-a)/n)*i;
        // std::cout<<h_x[i]<<std::endl;
    }

    // std::cout<<'1'<<std::endl;

    cudaMemcpy(x, h_x, n*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(sum, h_sum, n*sizeof(float), cudaMemcpyHostToDevice);

    // std::cout<<'1'<<std::endl;

    saxpy<<<n,1>>>(n, a, b, x, dx, sum);

    // std::cout<<'1'<<std::endl;
    // 'saxpy'의 함수를 1xn개의 gpu thread가 실행한다는 의미
    // saxpy<<<b, n>>>();
    // 에서 saxpy는 함수의 이름, b는 함수를 수행할 block의 개수ㅏ n은 다시 하나의 thread block 안에 몇 개의 thread가 존재하는지를 정하는 것.

    cudaMemcpy(h_sum, sum, n*sizeof(float), cudaMemcpyDeviceToHost);

    total = 0.;
    for(int i = 0; i<n; i++){
        float k = h_sum[i];
        total += k;
    }

    cudaFree(x);
    cudaFree(sum);

    free(h_x);
    free(h_sum);

    std::cout<<total<<std::endl;

    return 0;
}