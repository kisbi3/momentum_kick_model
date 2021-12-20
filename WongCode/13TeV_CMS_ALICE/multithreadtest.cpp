#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <fstream>
#include <thread>


void func1(double arr[]){
    for (int i = 0; i < 5; i++) {
        arr[i] += i;
        // std::cout<<arr[i];
    }
}

int main(){
    double arr1[5] = {0.,};
    double arr2[5];

    for (int i = 0; i < 5; i++) {
        arr2[i] = 1.;
    }


    for (int i = 0; i < 5; i++) {
        std::cout<<arr1[i]<<std::setw(10)<<arr2[i]<<std::endl;
    }
    std::thread t1(func1,arr1);
    std::thread t2(func1,arr2);
    t1.join();
    t2.join();
    // std::cout<<arr;
    std::cout<<std::endl;
    for (int i = 0; i < 5; i++) {
        std::cout<<arr1[i]<<std::setw(10)<<arr2[i]<<std::endl;
    }
    // t1.join();
    return 0;
}

// double func1(double num){
//     return num*2;
// }

// int main(){
//     double num, number;

//     num = 2.;

//     std::thread t1(func1,num);
//     t1.join();
//     std::cout<<num<<std::endl;
//     return 0;
// }