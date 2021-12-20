#include "cuda_function.cuh"
#include <cmath>


__device__ double fallof = .5;    //fall off parameter
__device__ double T = .5;    //Temperature, GeV
__device__ double q = .87;    //GeV
__device__ double m = 0.13957018;  //m == mpi
__device__ double mb = 0.13957018; //mb==mpi, GeV
__device__ double md = 1.;   //GeV
__device__ double sqrSnn = 200.;
__device__ double mp = 0.938272046; //Proton mass, GeV


//Jet Parameters
__device__ double Njet=10.;
__device__ double fj=1.;
__device__ double Tjet=0.7;
__device__ double sigmaphizero = .05;
__device__ double ma = 100.;

__device__ double Aridge;

__device__ double main_function(double x, double y){
    double dist = integralAridge(x, y);
    // printf("%f\n", dist);
    return dist;
    // return exp(x)+sqrt(y);
}

__device__ double integralNjet(double pt, double eta, double phi){
    double sigmaphi;
    double constant = Njet/(Tjet*(m+Tjet)*2*M_PI);
    sigmaphi = (sigmaphizero*sigmaphizero*ma*ma)/(ma*ma+pt*pt);
    return (constant/sigmaphi)*exp(((m-sqrt(m*m+pt*pt))/Tjet)-((phi*phi+eta*eta)/(2*sigmaphi)));
}

__device__ double lightcone(double pti, double yi){
    // printf("%f\n", sqrSnn);
    double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mp

    double squareroot=sqrt(m*m+pti*pti);
    // double yiabs = std::fabs(yi);
    return (squareroot/m)*exp(fabs(yi)-yb);
    // return exp(yiabs-yb);
}

__device__ double integralAridge(double pti, double yi){
    double x = lightcone(pti, yi);
    double squareroot=sqrt(m*m+pti*pti);
    // printf("%f\n", x);
    if(x>=1.){
        return 0.;
    }
    else{
        double fall = pow(1-x,0.5);
        printf("%f\n", fall);
        // return pti*pow(1-x,0.5)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
        return pti*sqrt(1.-x)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
        // return pti*exp(x)+yi;
        // return pti*pow(1.-x,a)+yi;
        // double power = pow(1.-x,6.);
        // return sqrt(1.-x);
    }
    
}

// __device__ double lightcone(double pti, double yi){    
//     double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mp

//     double squareroot=sqrt(m*m+pti*pti);
//     double yiabs = std::fabs(yi);
//     return (squareroot/m)*exp(yiabs-yb);
// }

__device__ double RidgeDisi(double pti, double yi){
    double x = lightcone(pti, yi);
    double squareroot=sqrt(m*m+pti*pti);


        if(x>=1.){
            return 0.;
        }
        else{
            return pow(1-x,fallof)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
        }
    
}


__device__ double RidgeDisf(double ptf, double etaf, double phif){
    double etajet = 0.;
    double ptisq = ptf*ptf-2*ptf*q*cos(phif)/cosh(etajet)+q*q/(cosh(etajet)*cosh(etajet));
    double pti;
    if(ptisq<0.0000000001){
        pti = 0.;
    }
    else{
        pti = sqrt(ptisq);
    }
    double E = sqrt(ptf*ptf*cosh(etaf)*cosh(etaf)+m*m);
    double Ei = sqrt(pti*pti+ptf*ptf*sinh(etaf)*sinh(etaf)+m*m);

    double yi = log((Ei+ptf*sinh(etaf))/(Ei-ptf*sinh(etaf)))/2;
    double yf = log((E+ptf*sinh(etaf))/(E-ptf*sinh(etaf)))/2;

    double x = lightcone(pti, yi);

    if (x>=1.){
        return 0.;
    }
    
    else{              
        return (Aridge*RidgeDisi(pti, yi))*sqrt(1.-((mb*mb)/((mb*mb+ptf*ptf)*cosh(yf)*cosh(yf))))*(E/Ei);
    }
    
}