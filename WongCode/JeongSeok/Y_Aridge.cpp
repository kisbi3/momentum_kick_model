#include <iostream>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#define PI 3.1415926535897932385
#define mb 0.938272046 // proton mass
#define m 0.13957018 // pion mass
using namespace std;

double Ridge(double pt, double y);
double Iy(double pt, double y0, double yn, int n);
double Ipt(double pt0, double ptn, double y0, double yn, int n);

int main()
{
    int n=100;
    double y0, yn, pt0, ptn;
    y0=-3.9;
    yn=3.9;
    pt0=1.15;
    ptn=10.;

    double Integral, Aridge;
    Integral=Ipt(pt0,ptn,y0,yn,n);
    Aridge=1./Integral;

    // cout << Aridge << endl;
    // cout << Aridge << endl;
    cout<<Integral<<setw(10)<<Aridge<<endl;
}

double Ridge(double pt, double y)
{
    double T, md, a, q, sqrt_sNN, yb; // constant
    T=0.50;
    md=1.0;
    a=0.5;
    q=1.0;
    sqrt_sNN=200.;
    yb=acosh(sqrt_sNN/(2.*m));
    
    double pti, x; // parameter
    pti=pt-q;
    x=sqrt(m*m+pti*pti)*exp(abs(y)-yb)/m; // lightcone variable

    double R_eta; // Ridge yi에 대한 식을 eta로 바꾼 것
    
    if (1-x<0.00001) {
        // cout<<'1'<<endl;
        R_eta=0.;
        // cout<<pti<<setw(10)<<y<<setw(10)<<R_eta<<endl;
        
    }
    else R_eta=pti*pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti); // Aridge을 포함하지 않은 eta에 대한 Ridge
    // cout<<pti<<setw(10)<<y<<setw(10)<<R_eta<<endl;
    return (R_eta);
}

// Aridge을 구하기 위한 적분(IAy->IApt)
double Iy(double pt, double y0, double yn, int n)
{
    int i;
    double y=0.,sum_mid=0.,dy;
    dy=(yn-y0)/n;

    for(i=1;i<n;i++)
    {
        y=y0+i*dy;
        sum_mid+=Ridge(pt,y);
    }

    double Integral;
    Integral=(Ridge(pt,y0)+Ridge(pt,yn)+2.*sum_mid)*dy/2.;

    return (Integral);
}

double Ipt(double pt0, double ptn, double y0, double yn, int n)
{
    int i,j;
    double pt=0.,sum_mid=0.,dpt, check;
    dpt=(ptn-pt0)/n;
    
    for(i=1;i<n;i++)
    {
        pt=pt0+i*dpt;
        //lightcone 계산 -> break
        check = Iy(pt,y0,yn,n);
        sum_mid+=check;
    }

    //ptn = pt-dpt;



    double Integral;
    Integral=2.*PI*(Iy(pt0,y0,yn,n)+Iy(ptn,y0,yn,n)+2.*sum_mid)*dpt/2.;

    return(Integral);
}
