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
    int n=1000;
    double y0, yn, pt0, ptn;
    y0=0.;
    yn=10.;
    pt0=0.;
    ptn=10.;

    double Integral, Aridge;
    Integral=Ipt(pt0,ptn,y0,yn,n);
    Aridge=1./Integral;

    // cout << Integral << endl;
    cout << Integral<<setw(15)<<Aridge << endl;
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

    double R_eta; // Ridge yi에 대한 식을 eta로 바꾼 것
    if (x>=1) R_eta=0.;
    else R_eta=2.*PI*pt*pow(1-x,a)*exp(-sqrt(m*m+pt*pt)/T)/sqrt(md*md+pt*pt);

    return (R_eta);
}

double lightcone(double pt, double y)
{
    double sqrt_sNN, yb; // constant
    sqrt_sNN=200.;
    yb=acosh(sqrt_sNN/(2.*mb));
    
    double x; // parameter
    x=sqrt(m*m+pt*pt)*exp(fabs(y)-yb)/m; // lightcone variable

    return (x);
}

// Aridge을 구하기 위한 적분(IAy->IApt)
double Iy(double pt, double y0, double yn, int n)
{
    int i;
    double y=0., sum_mid=0., dy, mid = 0.;
    dy=(yn-y0)/n;

    for(i=1;i<n;i++)
    {
        y=y0+i*dy;
        if (lightcone(pt,y)>1.)
        {
            sum_mid-=Ridge(pt,y-dy);
            // cout<<'1'<<endl;
            break;
        }
        else {
            // mid = Ridge(pt,y);
            // sum_mid+= mid;
            sum_mid+= Ridge(pt,y);
            // cout<<pt<<setw(10)<<y<<setw(15)<<lightcone(pt,y)<<setw(15)<<mid/(2.*PI)<<endl;
        }
    }

    // double y1=0., y2=0.;
    // for(i=1;i<n;i++)
    // {
    //     y1=y2+i*dy;
    //     if (lightcone(pt,y1)>1.) break;
    // }

    // cout << pt << setw(10) << y1 << setw(15) << lightcone(pt,y1-0.01) << endl;

    //lightcone이 y가 최대로 갔어도 1보다 크지 않을 경우 문제가 있음.
    double Integral;
    Integral=2.*(Ridge(pt,y0)+Ridge(pt,y-dy)+2.*sum_mid)*dy/2.;
    return (Integral);
}

double Ipt(double pt0, double ptn, double y0, double yn, int n)
{
    int i,j;
    double pt=0.,sum_mid=0.,dpt, checksum;
    dpt=(ptn-pt0)/n;
    
    for(i=1;i<n;i++)
    {
        pt=pt0+i*dpt;
        checksum =Iy(pt,y0,yn,n); 
        sum_mid+=checksum;
        // cout<<pt<<setw(15)<<Iy(pt,y0,yn,n)<<setw(15)<<sum_mid<<endl;
    }
    
    double Integral;
    Integral=(Iy(pt0,y0,yn,n)+Iy(ptn,y0,yn,n)+2.*sum_mid)*dpt/2.;;
    // cout<<pt<<setw(15)<<Iy(ptn,y0,yn,n)<<setw(15)<<sum_mid<<endl;
    return (Integral);
}
