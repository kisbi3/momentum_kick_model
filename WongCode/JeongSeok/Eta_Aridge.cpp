#include <iostream>
#include <stdio.h>
#include <math.h>
#define PI 3.1415926535897932385
#define mb 0.938272046 // proton mass
#define m 0.13957018 // pion mass
using namespace std;

double Ridge(double pt, double eta);
double Ieta(double pt, double eta0, double etan, int n);
double Ipt(double pt0, double ptn, double eta0, double etan, int n);

int main()
{
    int n=1000;
    double eta0, etan, pt0, ptn;
    eta0=-3.9;
    etan=3.9;
    pt0=0.15;
    ptn=4.;

    double Integral, Aridge;
    Integral=Ipt(pt0,ptn,eta0,etan,n);
    Aridge=1./Integral;

    cout << Aridge << endl;
}

double Ridge(double pt, double eta)
{
    double T, md, a, q, sqrt_sNN, yb; // constant
    T=0.50;
    md=1.0;
    a=0.5;
    q=1.0;
    sqrt_sNN=200.;
    yb=acosh(sqrt_sNN/(2.*mb));
    
    double pti, y_1, y_2, y, yi_1, yi_2, yi, x, Jacobian; // parameter
    pti=pt-q;
    y_1=sqrt(pow(pt*cosh(eta),2)+m*m);
    y_2=pt*sinh(eta);
    y=0.5*log((y_1+y_2)/(y_1-y_2));
    
    yi_1=sqrt(pow(pti*cosh(eta),2)+m*m);
    yi_2=pti*sinh(eta);
    yi=0.5*log((yi_1+yi_2)/(yi_1-yi_2));

    x=sqrt(m*m+pti*pti)*exp(abs(yi)-yb)/m; // lightcone variable

    Jacobian=sqrt(1-(m*m)/((m*m+pt*pt)*pow(cosh(y),2)));

    double R_eta; // Ridge yi에 대한 식을 eta로 바꾼 것
    
    if (x>=1.) R_eta=0.;
    else R_eta=Jacobian*pti*pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti); // Aridge을 포함하지 않은 eta에 대한 Ridge

    return (R_eta);
}

// Aridge을 구하기 위한 적분(IAy->IApt)
double Ieta(double pt, double eta0, double etan, int n)
{
    int i;
    double eta=0.,sum_mid=0.,deta;
    deta=(etan-eta0)/n;

    for(i=1;i<n;i++)
    {
        eta=eta0+i*deta;
        sum_mid+=Ridge(pt,eta);
    }

    double Integral;
    Integral=(Ridge(pt,eta0)+Ridge(pt,etan)+2.*sum_mid)*deta/2.;

    return (Integral);
}

double Ipt(double pt0, double ptn, double eta0, double etan, int n)
{
    int i,j;
    double pt=0.,sum_mid=0.,dpt;
    dpt=(ptn-pt0)/n;
    
    for(i=1;i<n;i++)
    {
        pt=pt0+i*dpt;
        sum_mid+=Ieta(pt,eta0,etan,n);
    }

    double Integral;
    Integral=2.*PI*(Ieta(pt0,eta0,etan,n)+Ieta(ptn,eta0,etan,n)+2.*sum_mid)*dpt/2.;

    return(Integral);
}