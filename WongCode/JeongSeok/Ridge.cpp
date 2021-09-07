#include <iostream>
#include <stdio.h>
#include <math.h>
#define PI 3.1415926535897932385
#define mb 0.938272046 // proton mass
#define m 0.13957018 // pion mass
using namespace std;

double Ridge(double pt, double eta);
double Ieta(double pt, double eta0, double etan, int n);

int main()
{
    FILE *fpt;

    fpt = fopen("Ridge.csv", "w+");
    fprintf(fpt,"pt, N_Rdige_dis\n");

    int n=1000;
    double pt, eta0, etan, N_Ridge_dis;
    eta0=-1.4;
    etan=1.4;

    if (pt=1) N_Ridge_dis=0.;
    else
    {
        for(pt=0.15;pt<=4;pt=pt+0.01)
        {
            N_Ridge_dis=Ieta(pt,eta0,etan,n);
            fprintf(fpt,"%f, %f\n",pt,N_Ridge_dis);
        }
        fclose(fpt);
    }
}

double Ridge(double pt, double eta)
{
    double T, md, a, sqrt_sNN, yb, Aridge; // constant
    T=0.5;
    md=1.0;
    a=0.5;
    sqrt_sNN=200.;
    yb=acosh(sqrt_sNN/(2.*mb));
    Aridge=0.1420836608;
    
    double pti, q;
    q=1.0;
    pti=pt-q;

    double y_1, y_2, y, E; // parameter
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

    double R_eta; // Ridge yi에 대한 식을 eta로 바꾼 것
    if (x>=1.) R_eta=0.;
    else R_eta=2.*PI*Aridge*Jacobian*(E/Ei)*pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti); // Aridge을 포함하지 않은 eta에 대한 Ridge

    return (R_eta);
}


double Ieta(double pt, double eta0, double etan, int n)
{
    double eta=0., sum_mid=0.;
    double deta=(etan-eta0)/n;
    int i;

    for(i=1;i<=n-1;i++)
    {
        eta=eta0+i*deta;
        sum_mid+=Ridge(pt,eta);
    }

    double Integral;
    Integral=(Ridge(pt,eta0)+Ridge(pt,etan)+2.*sum_mid)*deta/2.;

    return (Integral);
}

// double Ieta(double pt, double eta0, double etan, int n)
// {
//     int i;
//     double eta=0., mid4=0., mid2=0.;
//     double deta=(etan-eta0)/n;

//     for(i=1;i<=n-1;i+=2)
//     {
//         eta=eta0+i*deta;
//         mid4+=Ridge(pt,eta);
//     }
//     for(i=2;i<=n-2;i+=2)
//     {
//         eta=eta0+i*deta;
//         mid2+=Ridge(pt,eta);
//     }
//     return((Ridge(pt,eta0)+Ridge(pt,etan)+4.*mid4+2.*mid2)*deta/3.);
// }