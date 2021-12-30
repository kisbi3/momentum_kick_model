#include <stdio.h>
#include <iostream>
#include <time.h>
#include <thread>
#include <iomanip>


__device__ double lightcone(double pti, double yi, double sqrSnn, double mp, double m){
    // printf("%f\t%f\n", pti, yi);
    // printf("2");

    // double yi = rapidityinit(pti);
    double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mpi
    // double yb = __device__ double acosh(sqrSnn/(2.*mp));
    // printf("3\n");

    double squareroot=sqrt(m*m+pti*pti);
    // printf("%f\t%f\n", pti, yi);
    // printf("4\n");
    // double yiabs = abs(yi);

    // printf("5\n");
    // printf("%f, %f\n", pti, yi);
    double result = (squareroot/m)*exp(yi-yb);
    // printf("6\n");


    // printf("%f\t%f\n", yi, yiabs);
    // printf("%f, %f, %f\n", pti, yi, result);

    // std::cout<<pti<<std::setw(8)<<yi<<std::setw(15)<<exp(yiabs-yb)<<std::setw(15)<<squareroot/mb<<std::endl;
    // return (squareroot/m)*exp(yi-yb);
    // return (squareroot/m)*exp(fabs(yi)-yb);
    return result;

}

__device__ double integralAridge(double pti, double yi){
    // printf("%f\t%f\n", pti, yi);
    double sqrSnn = 200.;
    double mp = 0.938272046;
    double a = 0.5;
    double T = 0.5;
    double md = 1.;
    double m = 0.13957018;

    // printf("%f\t%f\n", pti, yi);
    // printf("1\n");
    double x = lightcone(pti, yi, sqrSnn, mp, m);
    // printf("%f, %f, %f\n", pti, yi, x);

    // double yi = rapidityinit(pti);
    // double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mpi
    // double squareroot=sqrt(m*m+pti*pti);
    // double yiabs = std::fabs(yi);
    // // std::cout<<pti<<std::setw(8)<<yi<<std::setw(15)<<exp(yiabs-yb)<<std::setw(15)<<squareroot/mb<<std::endl;
    // double x = (squareroot/m)*exp(yiabs-yb);


    // double squareroot=sqrt(m*m+pti*pti);
    // printf("%f\t%f\n", pti, yi);
    if(x>=1.){
        return 0.;
        // cal[index] = 0.;
    }
    else{
        return pti*pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
        // cal[index] = pti*pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
    }
}

// cal의 x를 pti, y를 yi로 놓고 계산하자. -> x : block / y : thread
__global__ void integrate(int n, double pti_start, double pti_end, double *pti, double yi_start, double yi_end, double *yi, double *sum){
    // printf("%d\n",blockIdx.x);
    // printf("%f\t%f\n", yi_end, yi_start);
    double dyi = ((yi_end-yi_start)/gridDim.x)*100;
    double dpti = (pti_end-pti_start)/blockDim.x;
    pti[blockIdx.x] = dpti*blockIdx.x+pti_start;
    yi[threadIdx.x] = dyi*threadIdx.x+yi_start;
    // pti = dpti*blockIdx.x+pti_start;
    // yi = dyi*threadIdx.x+yi_start;
    // printf("%d\n",blockIdx.x);
    // printf("%f\t%f\n", pti, yi);
    sum[blockIdx.x*blockDim.x+threadIdx.x] = integralAridge(pti[blockIdx.x], yi[threadIdx.x])*dyi*dpti;
    // printf("%d\n",blockIdx.x);
    // int index = blockIdx.x*blockDim.x+threadIdx.x;
    // cal[index] = integralAridge(cal[blockIdx.x], cal[blockIdx.y])
    // sum[index] = cal[index]*dx*dy;
}

// __global__ void integrate(double *pti, double dpti, double *yi, double dyi, double *sum){
//     // sum[blockIdx.x*gridDim.x+threadIdx.x] = integralAridge(pti[blockIdx.x], yi[threadIdx.x])*dyi*dpti;
//     sum[blockIdx.x+threadIdx.x*gridDim.x] = integralAridge(pti[blockIdx.x], yi[threadIdx.x])*dyi*dpti;
//     // printf("%f\t%f\n", pti[blockIdx.x], yi[threadIdx.x]);
// }

int main()
{
    clock_t start, finish;
    double duration;
    start = clock();
    using std::cout;
    using std::endl;
    using std::setw;

    // double dyi, dphii, sum, totalsum, phii, yi, dpti, pti, sum2, resultsum;
    double pti_start, pti_end, yi_start, yi_end, totalsum, Aridge;
    // double dyi, dpti;
    // int i, j, k, nyi, npti, nphii, check2;
    int nyi, npti;

    nyi = 100000;
    npti = 100000;

    pti_start = 0.;
    pti_end = 10.;
    yi_start = 0.;
    yi_end = 10.;

    // dyi = double((0.0+10.)/nyi);
    // dpti = double((0.+10.)/npti);

    // dyi = double ((yi_end-yi_start)/nyi);
    // dpti = double ((pti_end-pti_start)/npti);

    // dphii = double (M_PI+M_PI)/nphii;

    double *pti, *yi, *sum;
    // double *sum;

    // double *yi_cpu = (double*)malloc((nyi/100)*sizeof(double));
    // double *pti_cpu = (double*)malloc(npti*sizeof(double));
    double *sum_cpu = (double*)malloc(npti*(nyi/100)*sizeof(double));

    cudaMalloc((void**) &pti, npti*sizeof(double));
    cudaMalloc((void**) &yi, (nyi/100)*sizeof(double));
    cudaMalloc((void**) &sum, npti*(nyi/100)*sizeof(double));

    // double sum_cpu = (double*)malloc(npti*nyi*sizeof(double));


    // sum = 0.;
    // pti = 0.0;  //적분을 pt=0부터 하는것이 옳은가? 원통좌표계에서의 적분인데?
    // yi = 0.;    //0~4 적분한 후 x2할 것.

    // for(int i = 0;i<npti;i++){
    //     pti_cpu[i] = ((10.-0.)/npti)*i+pti_start;
    // }
    // cudaMemcpy(pti,pti_cpu, npti*sizeof(double), cudaMemcpyHostToDevice);

    totalsum = 0.;
    yi_end = 0.1;
    for(int k=0;k<100;k++){
        // for(int j=0;j<int(nyi/100);j++){
        //     yi_cpu[j] = ((10.-0.)/nyi)*j+yi_start;
        // }
        // cudaMemcpy(yi, yi_cpu, (nyi/100)*sizeof(double), cudaMemcpyHostToDevice);
        // cout<<"1"<<endl;
        // integrate<<<npti,int(nyi/100)>>>(pti, (pti_end-pti_start)/npti, yi, (yi_end-yi_start)/nyi, sum);
        printf("%f\t%f\n", yi_end, yi_start);
        integrate<<<npti,int(nyi/100)>>>(k, pti_start, pti_end, pti, yi_start, yi_end, yi, sum);
        printf("%f\t%f\n", yi_end, yi_start);
        // cout<<"222222"<<endl;
        cudaMemcpy(sum_cpu, sum, npti*(nyi/100)*sizeof(double), cudaMemcpyDeviceToHost);
        // cout<<"1"<<endl;
        for(int i=0; i<npti*(nyi/100);i++){
            totalsum += sum_cpu[i];
        }   //이 루프가 가장 느린듯
        printf("%f\t%f\t%f\n\n", yi_end, yi_start, totalsum);
        // cout<<"2"<<endl;
        yi_start += 0.1;
        yi_end += 0.1;
    }

    cudaFree(pti);
    cudaFree(yi);
    cudaFree(sum);

    // free(pti_cpu);
    // free(yi_cpu);
    free(sum_cpu);

    totalsum *= 4.*M_PI;
    Aridge = 1/totalsum;
    cout<<totalsum<<setw(20)<<Aridge<<endl;
    finish = clock();
    duration = (double)(finish-start)/CLOCKS_PER_SEC;
    cout<<"time : "<<duration<<" sec"<<endl;

    return 0;
}




