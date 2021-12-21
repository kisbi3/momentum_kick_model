#include <stdio.h>
//#include <cmath>
#ifndef cuda_function_cuh
#define cuda_function_cuh

// x : pti, y : yi
__device__ double main_function2(int function_number, double x, double y);
// x : ptf, y : etaf, z : phif
__device__ double main_function3(int function_number, double x, double y, double z);

//Jet
__device__ double integralNjet(double pt, double eta, double phi);

//lightcone
__device__ double lightcone(double pti, double yi);

//Aridge 적분
__device__ double integralAridge(double pti, double yi);

//Ridge distribution 안의 Aridge파트 적분
__device__ double RidgeDisi(double pti, double yi);

//Ridge distribution 적분
__device__ double RidgeDisf(double ptf, double etaf, double phif);



#endif