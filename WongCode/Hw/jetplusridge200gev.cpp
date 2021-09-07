#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<curses.h>
#include<iostream>
#define PI 3.1415926535897932385
#define mpi 0.13957018
#define mp 0.938272046

double jet(double pT, double eta, double phi);
double jetint(double pT, double etain, double etaout, double phi, int n);
double jetintint(double pT, double etain, double etaout, double phiin, double phiout, int n);

double ridge(double pT, double eta);
double ridgeint(double pT, double etain, double etaout, int n);

int main(){

    FILE *fpt;

    fpt = fopen("jetplusridge1.csv", "w+");
    fprintf(fpt,"pT, dis\n");

    int n=1000;
    double integral;
    double pT, eta, etain=-1.4, etaout=1.4,phiin=-1.,phiout=1.;

    for (pT=0;pT<=4.1;pT+=0.1){
        integral = (0.632*jetintint(pT, etain, etaout, phiin, phiout, n))+(/*2*PI**/2.53*ridgeint(pT, etain, etaout, n));
        fprintf(fpt,"%f, %f\n", pT, integral);
    }
    fclose(fpt);
}

// int main(){
//     int n=1000;
//     double integral;
//     double pT, eta, etain=-1.4, etaout=1.4,phiin=-1.,phiout=1.;

//     for (pT=0;pT<=4.1;pT+=0.1){
//         integral = (0.632*jetintint(pT, etain, etaout, phiin, phiout, n))+(/*2*PI**/2.53*ridgeint(pT, etain, etaout, n));
//         printf("{%f,%f},",pT,integral);
//     }
// }

double jet(double pT, double eta, double phi){
    double m=pow(mpi,2);
    double sigmaphi0=0.5;
    double ma=1.1;
    double Tjet=0.55;
    double Njet=0.75;
    double sigmaphi=sigmaphi0*ma/sqrt(pow(ma,2)+pow(pT,2));
    double jetcomponent;
    jetcomponent=Njet*((exp((mpi-sqrt(pow(m,2)+pow(pT,2))/Tjet))/(Tjet*(mpi+Tjet)))/(2*PI*pow(sigmaphi,2)))*exp(-(pow(phi,2)+pow(eta,2))/(2*pow(sigmaphi,2)));
    return(jetcomponent);
}

double ridge(double pT, double eta){
    double m=pow(mpi,2);
    double sqrtsnn=200.;
    double T=0.7;
    double a=0.5;
    double q=1.;
    double mmed=0.01;
    double md=1.0;
    
    double Aridge=0.023309;

    double pTinit=pT-q;
    double mT=sqrt(m+pow(pTinit,2));
    double mTf=sqrt(m+pow(pT,2));
    
    double xcomp=sqrt(m+pow((pTinit*cosh(eta)),2));
    double zcomp=pTinit*sinh(eta);
    double y=0.5*log((xcomp+zcomp)/(xcomp-zcomp));
    
    double xcompf=sqrt(m+pow((pT*cosh(eta)),2));
    double zcompf=pT*sinh(eta);
    double yf=0.5*log((xcompf+zcompf)/(xcompf-zcompf));
    
    double ybeam=acosh(sqrtsnn/(2*mp));
    double E=mTf*cosh(yf);
    double Ei=mT*cosh(y);

    double lightcone=exp(abs(y)-ybeam)*(mT/mpi);
    double jacobian=sqrt(1-(pow(mpi,2)/(pow(mTf,2)*pow(cosh(yf),2))));
    double fdist=pow((1-lightcone),a)*(exp(-mT/T)/sqrt(pow(md,2)+pow(pTinit,2)))*(E/Ei);
    double ridgecomponent;

    if (lightcone<=1){
        ridgecomponent=Aridge*jacobian*fdist;
    }
    else{
        ridgecomponent=0;
    }

    return(ridgecomponent);
}

double jetint(double pT, double phi, double etain, double etaout, int n){
    int i;
    double eta=0, sum_mid=0., h;

    h=(etaout-etain)/n;

    for(i=1;i<=n-1;i++){
        eta = etain+i*h;
        sum_mid += jet(pT,eta,phi);
    }
    return((jet(etain,pT,phi)+jet(etaout,pT,phi)+2.*sum_mid)*h/2.);
}

double jetintint(double pT, double etain, double etaout, double phiin, double phiout, int n){
    int i;
    double phi=0, sum_mid=0., phih;

    phih=(phiout-phiin)/n;

    for(i=1;i<=n-1;i++){
        phi = phiin+i*phih;
        sum_mid += jetint(pT,phi,etain,etaout,n);
    }
    return((jetint(phiin,pT,etain,etaout,n)+jetint(phiout,pT,etain,etaout,n)+2.*sum_mid)*phih/2.);
}

double ridgeint(double pT, double etain, double etaout, int n){
    int i;
    double eta=0, sum_mid=0., etah;

    etah=(etaout-etain)/n;

    for(i=1;i<=n-1;i++){
        eta = etain+i*etah;
        sum_mid += ridge(pT,eta);
    }
    return((ridge(etain,pT)+ridge(etaout,pT)+2.*sum_mid)*etah/2.);
}