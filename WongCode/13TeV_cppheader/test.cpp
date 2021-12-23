#include <iostream>
#include <cmath>
#include <iomanip>
#include "integral.hpp"
#include "function.hpp"
// #include "headerfiles/function.hpp"


// double func1(double x, double y, double z){
//     return x*x+y*y+z*z;
// }

void func1(int num){
    Aridge Ar;
    Ar.integralAridge();
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
    // Aridge A;
    // cuda_integral cuda;
    // double result = i.secondintegral(func1, 0., 1., 100, 0., 1., 100);

    // cout<<result<<endl;
    // cout << "CPU : "<< integrate.secondintegral(func1, 0., 1., 10000, 0., 1., 10000) << endl;
    // cout << "stdarg : "<< integrate.quadintegral(3, 0., 1., .5, 3., 4.4, .2, 6., 7., .5) << endl;
    cout << "stdarg : "<< integrate.quadintegral(func1, 3, 0., 3., .01, 0., 3., .01, 0., 3., .01) << endl;

    // cout<< testint(func1, 0., 1., 100, 0., 1., 100) << endl;

    return 0;
}