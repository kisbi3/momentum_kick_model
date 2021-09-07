#include <iostream>
#include <stdio.h>
#include <math.h>
#define PI 3.1415926535897932385
#define mb 0.938272046 // proton mass
#define m 0.13957018 // pion mass
using namespace std;

double Ridge(double pt, double y);

int main()
{
    FILE *fpt;

    double pt, y, Nor;
    pt=4.;

    fpt=fopen("4_Nor_ini_par_mom_dis.csv", "w+");
    fprintf(fpt,"y, dF/dyptdpt, x\n");

    for(y=-10.;y<=10.;y=y+0.01)
    {
        Nor=Ridge(pt,y);
        fprintf(fpt,"%f, %f\n", y, Nor);
    }
    fclose(fpt);
}

double Ridge(double pt, double y)
{
    double T, md, a, sqrt_sNN, yb; // constant
    T=0.50;
    md=1.0;
    a=0.5;
    sqrt_sNN=200.;
    yb=acosh(sqrt_sNN/(2.*mb));
    
    double x; // parameter
    x=sqrt(m*m+pt*pt)*exp(abs(y)-yb)/m; // lightcone variable

    double dF, Aridge; // Ridge yi에 대한 식을 eta로 바꾼 것
    Aridge=0.443086;

    if (x>=1.) dF=0.;
    else dF=2*PI*Aridge*pow(1-x,a)*exp(-sqrt(m*m+pt*pt)/T)/sqrt(md*md+pt*pt); // Aridge을 포함하지 않은 eta에 대한 Ridge

    return (dF);
}