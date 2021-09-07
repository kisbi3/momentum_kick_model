#include <iostream>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#define PI 3.1415926535897932385
#define mb 0.938272046 // proton mass
#define m 0.13957018 // pion mass
using namespace std;

double Ridge(double pt, double y);
double lightcone(double pt, double y);
double Iy(double pt, double y0, double yn, int n);
double Ipt(double pt0, double ptn, double y0, double yn, int n);

int main()
{
    int n=10000;
    double y0, yn, pt0, ptn;
    y0=0.; // -6에서 6까지한 거
    yn=6.;
    pt0=0.;
    ptn=10.;

    double Integral, Aridge;
    Integral=Ipt(pt0,ptn,y0,yn,n);
    Aridge=1./Integral;

    // cout << Integral << endl;
    cout << Aridge << endl;
}

double Ridge(double pt, double y)
{
    double T, md, a; // constant
    T=0.50;
    md=1.0;
    a=0.5;

    double sqrt_sNN, yb; // constant
    sqrt_sNN=200.;
    yb=acosh(sqrt_sNN/(2.*mb));

    double x;
    x=sqrt(m*m+pt*pt)*exp(abs(y)-yb)/m;

    double R_y; // Ridge y에 대한 식
    R_y=2.*PI*pt*pow(1-x,a)*exp(-sqrt(m*m+pt*pt)/T)/sqrt(md*md+pt*pt);

    return (R_y);
}

double lightcone(double pt, double y)
{
    double sqrt_sNN, yb; // constant
    sqrt_sNN=200.;
    yb=acosh(sqrt_sNN/(2.*mb));
    
    double x; // parameter
    x=sqrt(m*m+pt*pt)*exp(abs(y)-yb)/m; // lightcone variable

    return (x);
}

// Aridge을 구하기 위한 적분(IAy->IApt)
double Iy(double pt, double y0, double yn, int n)
{
    int i;
    double y=0., sum_mid=0., dy;
    dy=(yn-y0)/n;

    for(i=1;i<n;i++)
    {
        y=y0+i*dy;
        if (lightcone(pt,y)>1.)
        {
            sum_mid-=Ridge(pt,y-dy);
            break;
        }
        else sum_mid+=Ridge(pt,y);
    }

    // cout << pt << setw(10) << y1 << setw(15) << lightcone(pt,y1-0.01) << endl;

    double Integral;
    Integral=2.*(Ridge(pt,y0)+Ridge(pt,y-dy)+2.*sum_mid)*dy/2.;
    return (Integral);
}

double Ipt(double pt0, double ptn, double y0, double yn, int n)
{
    int i,j;
    double pt=0.,sum_mid=0.,dpt;
    dpt=(ptn-pt0)/n;
    
    for(i=1;i<n;i++)
    {
        pt=pt0+i*dpt;
        if (Iy(pt,y0,yn,n)<0.00001)
        {
            sum_mid-=Iy(pt-dpt,y0,yn,n);
            break;
        }
        else sum_mid+=Iy(pt,y0,yn,n);
    }

    // cout << Iy(pt-dpt,y0,yn,n) << endl;
    
    double Integral;
    Integral=(Iy(pt0,y0,yn,n)+Iy(pt-dpt,y0,yn,n)+2.*sum_mid)*dpt/2.;
    return (Integral);
}
