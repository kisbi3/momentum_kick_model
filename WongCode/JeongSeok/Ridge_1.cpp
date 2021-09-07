#include <iostream>
#include <stdio.h>
#include <math.h>
#define PI 3.1415926535897932385
#define mb 0.938272046 // proton mass
#define m 0.13957018 // pion mass
using namespace std;

double Ridge(double pt, double eta, double phi);
double lightcone(double pt, double eta, double phi);
double IReta(double pt, double phi, double R_eta0, double etan, int n);
double IRphi(double pt, double R_eta0, double etan, double phi0, double phin, int n);

int main()
{
    FILE *fpt;

    fpt = fopen("Ridge_1.csv", "w+");
    fprintf(fpt,"pt, N_Rdige_dis\n");

    int n=100;
    double pt, R_eta0, etan, phi0, phin, N_Ridge_dis;
    R_eta0=-1.4; // -1.4에서 1.4까지한거
    etan=1.4;
    phi0=-1.;
    phin=1.;

    for(pt=0.15;pt<=4;pt=pt+0.01)
        {
            N_Ridge_dis=IRphi(pt,R_eta0,etan,phi0,phin,n);
            fprintf(fpt,"%f, %f\n",pt,N_Ridge_dis);
        }
        fclose(fpt);
}

double Ridge(double pt, double eta, double phi)
{
    double T, md, a, sqrt_sNN, yb;
    T=0.5;
    md=1.0;
    a=0.5;
    sqrt_sNN=200.;
    yb=acosh(sqrt_sNN/(2.*mb));
    
    double pti, q;
    q=1.0;
    pti=sqrt(pt*pt-2.*pt*q*cos(phi)+q*q);

    double y_1, y_2, y, E;
    y_1=sqrt(pow(pt*cosh(eta),2)+m*m);
    y_2=pt*sinh(eta);
    y=0.5*log((y_1+y_2)/(y_1-y_2));
    E=sqrt(m*m+pt*pt)*cosh(y);
    
    double yi_1, yi_2, yi, Ei;
    yi_1=sqrt(pow(pti*cosh(eta),2)+m*m);
    yi_2=pti*sinh(eta);
    yi=0.5*log((yi_1+yi_2)/(yi_1-yi_2));
    Ei=sqrt(m*m+pti*pti)*cosh(yi);

    double x;
    x=sqrt(m*m+pti*pti)*exp(abs(yi)-yb)/m; // lightcone variable

    double Jacobian;
    Jacobian=sqrt(1-(m*m)/((m*m+pt*pt)*pow(cosh(y),2)));

    double Aridge;
    Aridge=0.142039;

    double R_eta; // Ridge yi에 대한 식을 eta로 바꾼 것
    R_eta=Aridge*(E/Ei)*Jacobian*pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti); // Aridge을 포함하지 않은 eta에 대한 Ridge
    //*(E/Ei)
    return (R_eta);
}

double lightcone(double pt, double eta, double phi)
{
    double sqrt_sNN, yb; // constant
    sqrt_sNN=200.;
    yb=acosh(sqrt_sNN/(2.*mb));

    double pti, q;
    q=1.0;
    pti=sqrt(pt*pt-2.*pt*q*cos(phi)+q*q);

    double yi_1, yi_2, yi;
    yi_1=sqrt(pow(pti*cosh(eta),2)+m*m);
    yi_2=pti*sinh(eta);
    yi=0.5*log((yi_1+yi_2)/(yi_1-yi_2));
    
    double x; // parameter
    x=sqrt(m*m+pti*pti)*exp(abs(yi)-yb)/m; // lightcone variable

    return (x);
}

double IReta(double pt, double phi, double R_eta0, double etan, int n)
{
    double eta=0., sum_mid=0.;
    double deta=(etan-R_eta0)/n;
    int i;

    for(i=1;i<n;i++)
    {
        eta=R_eta0+i*deta;
        if (lightcone(pt,eta,phi)>1.)
        {
            sum_mid-=Ridge(pt,eta-deta,phi);
            break;
        }
        else sum_mid+=Ridge(pt,eta,phi);
    }

    double Integral;
    Integral=(Ridge(pt,R_eta0,phi)+Ridge(pt,eta-deta,phi)+2.*sum_mid)*deta/2.;

    return (Integral);
}

double IRphi(double pt, double R_eta0, double etan, double phi0, double phin, int n)
{
    double phi=0., sum_mid=0.;
    double dphi=(phin-phi0)/n;
    int i;

    for(i=1;i<n;i++)
    {
        phi=phi0+i*dphi;
        if (IReta(pt,phi,R_eta0,etan,n)<0.00001)
        {
            sum_mid-=IReta(pt,phi-dphi,R_eta0,etan,n);
            break;
        }
        else sum_mid+=IReta(pt,phi,R_eta0,etan,n);
    }
    
    double Integral;
    Integral=(IReta(pt,phi0,R_eta0,etan,n)+IReta(pt,phi-dphi,R_eta0,etan,n)+2.*sum_mid)*dphi/2.;
    return (Integral);
}