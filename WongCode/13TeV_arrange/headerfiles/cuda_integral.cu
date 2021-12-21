// #include "integral.hpp"
#include "cuda_function.cuh"
#include "cuda_integral.cuh"

// static const int blockSize = 1024;

// __device__ double inside_cuda_function(double (*func)(double, double), double x, double y){
//     return x*x+y*y;
// }

__global__ void integrate(double x_st, double x_end, double *x, double y_st, double y_end, double *y, double *sum){
    double dx = (x_end-x_st)/gridDim.x;
    double dy = (y_end-y_st)/blockDim.x;

    sum[blockIdx.x*blockDim.x+threadIdx.x] = main_function(dx*blockIdx.x+x_st, dy*threadIdx.x+y_st)*dx*dy;
    // __syncthreads();
}

__global__ void sumCommMultiBlock(double *gArr, int arraySize, double *gOut) {
    static const int blockSize = 1024;
    
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x*blockDim.x;
    const int gridSize = blockDim.x*gridDim.x;
    double sum = 0;
    // printf("%f\n", gArr[threadIdx.x]);
    for (int i = gthIdx; i < arraySize; i += gridSize)
        sum += gArr[i];
    __shared__ double shArr[blockSize];
    shArr[thIdx] = sum;
    __syncthreads();
    for (int size = blockDim.x/2; size>0; size/=2) { //uniform
        if (thIdx<size)
            shArr[thIdx] += shArr[thIdx+size];
        __syncthreads();
    }
    if (thIdx == 0)
        gOut[blockIdx.x] = shArr[0];
}



extern "C" double cuda_twodimsumarray(double *gArr, int bin_x, int bin_y, double *gOut){
    double totaldist = 0.;
    sumCommMultiBlock<<<bin_x,1024>>>(gArr, bin_x*bin_y, gOut);
    sumCommMultiBlock<<<1,1024>>>(gOut, bin_x, gOut);
    cudaDeviceSynchronize();
    cudaMemcpy(&totaldist, gOut, sizeof(double), cudaMemcpyDeviceToHost);
    return totaldist;
}

// double cuda_function(double (*func)(double, double), double x, double y){

    
//     return func(x, y);
// }

extern "C" double cuda_secondintegral(double x_st, double x_end, int bin_x, double y_st, double y_end, int bin_y){

    double totalsum = 0.;
    double *x, *y, *sum, *sumdist;
    cudaMalloc((void**) &x, bin_x*sizeof(double));
    cudaMalloc((void**) &sumdist, bin_x*sizeof(double));

    // if(bin_y>1000){
    //     #define over
    // }
    
    // #define over (bin_y>1000)
    // // printf("\n%d", bin_y);
    // #ifdef over
    int thread_num = 800;
    if(bin_y>thread_num){
        #define overthread
    }
    // if(bin_y>thread_num){
    // #if bin_y > thread_num
    #ifdef overthread
        // printf("\n1");
        int iteration = bin_y/thread_num;
        double dy = (y_end-y_st)/iteration;
        y_end = dy;
        
        double *sum_cpu = (double*)malloc(bin_x*(bin_y/iteration)*sizeof(double));
        cudaMalloc((void**) &y, (bin_y/iteration)*sizeof(double));
        cudaMalloc((void**) &sum, bin_x*(bin_y/iteration)*sizeof(double));    
        for(int k=0;k<iteration;k++){
            // double totaldist = 0.;
            // integrate<<<bin_x,1000>>>(func, x_st, x_end, x, y_st, y_end, y, sum);
            integrate<<<bin_x,thread_num>>>(x_st, x_end, x, y_st, y_end, y, sum);
            // sumCommMultiBlock<<<bin_x,1024>>>(sum, bin_x*(bin_y/100), sumdist);
            // sumCommMultiBlock<<<1,1024>>>(sumdist, bin_x, sumdist);
            // cudaDeviceSynchronize();
            // cudaMemcpy(&totaldist, sumdist, sizeof(double), cudaMemcpyDeviceToHost);

            double totaldist = cuda_twodimsumarray(sum, bin_x, bin_y/iteration, sumdist);

            // printf("%f\n", totaldist);
            totalsum += totaldist;
            y_st += dy;
            y_end += dy;
        }
        free(sum_cpu);
    // }
    #else 
    // else{
        // printf("111");
        double *sum_cpu = (double*)malloc(bin_x*bin_y*sizeof(double));
        cudaMalloc((void**) &y, bin_y*sizeof(double));
        cudaMalloc((void**) &sum, bin_x*bin_y*sizeof(double)); 
        integrate<<<bin_x,bin_y>>>(x_st, x_end, x, y_st, y_end, y, sum);
        totalsum = cuda_twodimsumarray(sum, bin_x, bin_y, sumdist);
        free(sum_cpu);
    // }
    #endif



    //Error Message!!
    cudaError_t err = cudaGetLastError();

    if ( err != cudaSuccess )
    {
       printf("CUDA Error: %s\n", cudaGetErrorString(err));       

       // Possibly: exit(-1) if program cannot continue....
    }

    cudaFree(x);
    cudaFree(y);
    cudaFree(sum);
    cudaFree(sumdist);

    // free(sum_cpu);

    return totalsum;
}

// extern "C" double cuda_thirdintegral(double x_st, double x_end, int bin_x, double y_st, double y_end, int bin_y){

//     double totalsum = 0.;
//     double *x, *y, *z, *sum, *sumdist;
//     cudaMalloc((void**) &x, bin_x*sizeof(double));
//     cudaMalloc((void**) &sumdist, bin_x*sizeof(double));

//     // if(bin_y>1000){
//     //     #define over
//     // }
    
//     // #define over (bin_y>1000)
//     // // printf("\n%d", bin_y);
//     // #ifdef over
//     int thread_num = 800;
//     if(bin_y>thread_num){
//         #define overthread
//     }
//     // if(bin_y>thread_num){
//     // #if bin_y > thread_num
//     #ifdef overthread
//         // printf("\n1");
//         int iteration = bin_y/thread_num;
//         double dy = (y_end-y_st)/iteration;
//         y_end = dy;
        
//         double *sum_cpu = (double*)malloc(bin_x*(bin_y/iteration)*sizeof(double));
//         cudaMalloc((void**) &y, (bin_y/iteration)*sizeof(double));
//         cudaMalloc((void**) &sum, bin_x*(bin_y/iteration)*sizeof(double));    
//         for(int k=0;k<iteration;k++){
//             // double totaldist = 0.;
//             // integrate<<<bin_x,1000>>>(func, x_st, x_end, x, y_st, y_end, y, sum);
//             integrate<<<bin_x,thread_num>>>(x_st, x_end, x, y_st, y_end, y, sum);
//             // sumCommMultiBlock<<<bin_x,1024>>>(sum, bin_x*(bin_y/100), sumdist);
//             // sumCommMultiBlock<<<1,1024>>>(sumdist, bin_x, sumdist);
//             // cudaDeviceSynchronize();
//             // cudaMemcpy(&totaldist, sumdist, sizeof(double), cudaMemcpyDeviceToHost);

//             double totaldist = cuda_twodimsumarray(sum, bin_x, bin_y/iteration, sumdist);

//             // printf("%f\n", totaldist);
//             totalsum += totaldist;
//             y_st += dy;
//             y_end += dy;
//         }
//         free(sum_cpu);
//     // }
//     #else 
//     // else{
//         // printf("111");
//         double *sum_cpu = (double*)malloc(bin_x*bin_y*sizeof(double));
//         cudaMalloc((void**) &y, bin_y*sizeof(double));
//         cudaMalloc((void**) &sum, bin_x*bin_y*sizeof(double)); 
//         integrate<<<bin_x,bin_y>>>(x_st, x_end, x, y_st, y_end, y, sum);
//         totalsum = cuda_twodimsumarray(sum, bin_x, bin_y, sumdist);
//         free(sum_cpu);
//     // }
//     #endif



//     //Error Message!!
//     cudaError_t err = cudaGetLastError();

//     if ( err != cudaSuccess )
//     {
//        printf("CUDA Error: %s\n", cudaGetErrorString(err));       

//        // Possibly: exit(-1) if program cannot continue....
//     }

//     cudaFree(x);
//     cudaFree(y);
//     cudaFree(sum);
//     cudaFree(sumdist);

//     // free(sum_cpu);

//     return totalsum;
// }