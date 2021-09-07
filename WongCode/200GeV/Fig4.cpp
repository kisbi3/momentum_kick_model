#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstring>

double a = 0.5; //fall off parameter
double T = 0.5; //Temperature, Gev
double q = 1.;  //GeV
double fRN = 4.;
double md = 1.; //GeV
double mb = 0.13957018; //m==mb==mpi, GeV
double sqrSnn = 200; //GeV
double mp = 0.938; //Proton mass, GeV




double lightcone(double pti, double yi){
    // double yi = rapidityinit(pti);
    double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mpi
    double squareroot=sqrt(mb*mb+pti*pti);
    double yiabs = std::fabs(yi);
    // std::cout<<pti<<std::setw(8)<<yi<<std::setw(15)<<exp(yiabs-yb)<<std::setw(15)<<squareroot/mb<<std::endl;
    return (squareroot/mb)*exp(yiabs-yb);
}

//pt, y 대입할 함수, pt는 없어야 한다.
double integralAridge(double pti, double yi){
    // std::cout<<pti<<std::setw(7)<<yi<<std::setw(9)<<phii<<std::endl;       
    double x = lightcone(pti, yi);
    double squareroot=sqrt(mb*mb+pti*pti);
    double Ei = sqrt(mb*mb+(pti)*(pti))*cosh(yi);
    // std::cout<<pti<<std::setw(10)<<yi<<std::setw(10)<<phii<<std::setw(10)<<x<<std::endl;
    if (x>=1.){
        return 0.;
    }
    return 0.1420386444451514*2*M_PI*pow((1-x),a)*((exp(-squareroot/T))/(sqrt(md*md+pti*pti)));
    //
}


int main()
{
    using std::cout;
    using std::endl;
    using std::setw;
    double dyi, dpti, pti, yi;
    int i, k, n;

    double result1[10][1000] = {0.}, result2[10][1000] = {0.}, sum[10] = {0.};

    FILE *fpt;
    fpt = fopen("Test.csv", "w+");
    fprintf(fpt,"pt, y=0, y=1, y=2, y=3, y=4\n");
    


    n = 1000000;
    dpti = double ((4.-0.)/n);
    yi = 0.;

    pti = 0.;
    for(k=1;k<=n+1;k+=1){
        fprintf(fpt,"%f, %f, %f, %f, %f, %f\n", pti, integralAridge(pti, 0.), integralAridge(pti, 1.), integralAridge(pti, 2.), integralAridge(pti, 3.), integralAridge(pti, 4.));
        pti += dpti;
    }

    fpt = fopen("Test2.csv", "w+");
    fprintf(fpt,"yi, pt=0, pt=1, pt=2, pt=3, pt=4\n");
    


    n = 1000000;
    dyi = double ((6.+6.)/n);
    yi = -6.;
    pti = 0.1;

    yi = -6.;
    for(k=1;k<=n+1;k+=1){
        fprintf(fpt,"%f, %f, %f, %f, %f, %f\n", yi, integralAridge(0., yi), integralAridge(1., yi), integralAridge(2., yi), integralAridge(3., yi), integralAridge(4., yi));
        yi += dyi;
    }


    fclose(fpt);


    return 0;
}
