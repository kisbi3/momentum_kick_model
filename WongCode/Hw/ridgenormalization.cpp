#include<stdio.h>
#include<math.h>
#include<stdlib.h>
// #include<curses.h>
#include<iostream>
#define PI 3.1415926535897932385
#define mpi 0.13957018
#define mp 0.938272046

double ridge(double pT, double eta);
double ridgeint(double pT, double etain, double etaout, int n);
double ridgeintint(double pTin, double pTout, double etain, double etaout, int n);

int main(){
    int n=1000;
    double Aridge;
    double pT, eta, etain=-8., etaout=8., pTin=0., pTout=5.;
    Aridge = 1/(2*PI*ridgeintint(pTin, pTout, etain, etaout, n));
    printf("%f\n",Aridge);
}

double ridge(double pT, double eta){
    double m=pow(mpi,2);
    double sqrtsnn=200.;
    double T=0.5;
    double a=0.5;
    double q=1.;
    double mmed=0.01;
    double md=1.0;
    
    double pTinit=pT-q;
    double mT=sqrt(m+pow(pTinit,2));
    double mTf=sqrt(m+pow(pT,2));
    
    double xcomp=sqrt(m+pow((pTinit*cosh(eta)),2));
    double zcomp=pTinit*sinh(eta);
    double y=0.5*log((xcomp+zcomp)/(xcomp-zcomp));
    
    // double xcompf=sqrt(m+pow((pT*cosh(eta)),2));
    // double zcompf=pT*sinh(eta);
    // double yf=0.5*log((xcompf+zcompf)/(xcompf-zcompf));
    
    double ybeam=acosh(sqrtsnn/(2*mp));
    // double E=mTf*cosh(yf);
    // double Ei=mT*cosh(y);

    double lightcone=exp(abs(y)-ybeam)*(mT/mpi);
    // double jacobian=sqrt(1-(pow(mpi,2)/(pow(mTf,2)*pow(cosh(yf),2))));
    double fdist=pow((1-lightcone),a)*(exp(-mT/T)/sqrt(pow(md,2)+pow(pTinit,2)));
    double ridgecomponent;

    if (lightcone<=1){
        ridgecomponent=/*jacobian*/pT*fdist;
    }
    else{
        ridgecomponent=0;
    }

    return(ridgecomponent);
}

double ridgeint(double pT, double etain, double etaout, int n){
    int i;
    double eta=0, sum_mid=0., etah;

    etah=(etaout-etain)/n;

    for(i=1;i<=n-1;i++){
        eta = etain+i*etah;
        sum_mid += ridge(pT,eta);
    }
    return((ridge(pT,etain)+ridge(pT,etaout)+2.*sum_mid)*etah/2.);
}

double ridgeintint(double pTin, double pTout, double etain, double etaout, int n){
    int i;
    double pT=0, sum_mid=0., pTh;

    pTh=(pTout-pTin)/n;

    for(i=1;i<=n-1;i++){
        pT = pTin+i*pTh;
        sum_mid += ridgeint(pT,etain,etaout,n);
    }
    return((ridgeint(pTin,etain,etaout,n)+ridgeint(pTout,etain,etaout,n)+2.*sum_mid)*pTh/2.);
}