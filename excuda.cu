// #include <stdio.h>

// __global__
// void saxpy(int n, float a, float *x, float *y)
// {
//   int i = blockIdx.x*blockDim.x + threadIdx.x;
//   if (i < n) y[i] = a*x[i] + y[i];
// }

// int main(void)
// {
//   int N = 1<<20;
//   float *x, *y, *d_x, *d_y;
//   x = (float*)malloc(N*sizeof(float));
//   y = (float*)malloc(N*sizeof(float));

//   cudaMalloc(&d_x, N*sizeof(float)); 
//   cudaMalloc(&d_y, N*sizeof(float));

//   for (int i = 0; i < N; i++) {
//     x[i] = 1.0f;
//     y[i] = 2.0f;
//   }

//   cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
//   cudaMemcpy(d_y, y, N*sizeof(float), cudaMemcpyHostToDevice);

//   // Perform SAXPY on 1M elements
//   saxpy<<<(N+255)/256, 256>>>(N, 2.0f, d_x, d_y);

//   cudaMemcpy(y, d_y, N*sizeof(float), cudaMemcpyDeviceToHost);

//   float maxError = 0.0f;
//   for (int i = 0; i < N; i++)
//     maxError = max(maxError, abs(y[i]-4.0f));
//   printf("Max error: %f\n", maxError);

//   cudaFree(d_x);
//   cudaFree(d_y);
//   free(x);
//   free(y);
// }

#include <iostream>

__global__ void saxpy(int n, float a, float *__restrict__ x, float *__restrict__ y)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if( i < n )
    y[i] = a * x[i] + y[i];
}

int main()
{
  std::cout<<'1'<<std::endl;
  int N = 1 << 16;
  
  int size = N * sizeof(float);
  float *h_x = (float*)malloc(size);
  float *h_y = (float*)malloc(size);
  
  float *d_x;
  float *d_y;
  
  cudaMalloc((void**) &d_x, size);
  cudaMalloc((void**) &d_y, size);
  
  std::cout<<'1'<<std::endl;
  for(int i=0; i < N; i++)
  {
    h_x[i] = 2.0;
    h_y[i] = 2.0;
  }
  
  cudaMemcpy(d_x, h_x, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, h_y, size, cudaMemcpyHostToDevice);
  
  std::cout<<'1'<<std::endl;
  saxpy<<<256, 256>>>(N, 2.0, d_x, d_y);
  
  cudaMemcpy(h_y, d_y, size, cudaMemcpyDeviceToHost);

  std::cout<<*h_y<<std::endl;
  
  cudaFree(d_x);
  cudaFree(d_y);
  
  free(h_x);
  free(h_y);

  
  
  return 0;
}