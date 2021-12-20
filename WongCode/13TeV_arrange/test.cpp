#include <iostream>
#include <cmath>
#include <iomanip>
#include "headerfiles/integral.hpp"
#include "headerfiles/function.hpp"
#include "headerfiles/cuda_integral.cuh"


double func1(double x, double y){
    return x*x+y*y;
}

double cuda_function(double x, double y){
    return x*x+y*y;
}

// __host__ __device__ double cuda_func(double x, double y){
//     return x*x+y*y;
// }

//function들을 모두 cuda형식으로 바꿔서 진행해보자.

int main()
{
    using std::cout;
    using std::endl;
    integral integrate;
    Aridge A;
    // cuda_integral cuda;
    // double result = i.secondintegral(func1, 0., 1., 100, 0., 1., 100);

    // cout<<result<<endl;
    cout << "CPU : "<< integrate.secondintegral(func1, 0., 1., 10000, 0., 1., 10000) << endl;
    // cout << cuda.secondintegral(cuda_function, 0., 1., 100000, 0., 1., 100000) << endl;
    // cout << cuda_secondintegral(cuda_func, 0., 1., 100000, 0., 1., 100000) << endl;
    double totalsum = cuda_secondintegral(0., 10., 10000, 0., 10., 10000);
    totalsum *= 4.*M_PI;
    cout << "GPU : "<< totalsum << std::setw(20)<< 1/totalsum << endl;

    // cout<< testint(func1, 0., 1., 100, 0., 1., 100) << endl;

    return 0;
}