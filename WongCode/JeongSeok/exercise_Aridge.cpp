#include <iostream>
#include <stdio.h>
#include <math.h>
#define PI 3.1415926535897932385
#define mb 0.938272046 // proton mass
#define m 0.13957018 // pion mass
using namespace std;

double Ridge(double pt, double eta);

int main()
{
    // double pt, eta;
    // pt = 10.;
    // eta = 2.;

    cout << Ridge(1.,18.) << endl;

    // double T, md, a, q, sqrt_sNN, yb; // constant
    // T=0.65;
    // md=1.5;
    // a=0.5;
    // q=0.3;
    // sqrt_sNN=200.;
    // yb=acosh(sqrt_sNN/(2.*mb));
    
    // double pti, y_1, y_2, y, yi_1, yi_2, yi, x, Jacobian; // parameter
    // pti=pt-q;
    // y_1=sqrt(pow(pt*cosh(eta),2)+m*m);
    // y_2=pt*sinh(eta);
    // y=0.5*log((y_1+y_2)/(y_1-y_2));

    // // cout << y << endl;
    
    // yi_1=sqrt(pow(pti*cosh(eta),2)+m*m);
    // yi_2=pti*sinh(eta);
    // yi=0.5*log((yi_1+yi_2)/(yi_1-yi_2));

    // // cout << yi << endl;

    // x=sqrt(m*m+pti*pti)*exp(abs(yi)-yb)/mb; // lightcone variable

    // cout << x << endl;
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

    // cout << y << endl;
    
    yi_1=sqrt(pow(pti*cosh(eta),2)+m*m);
    yi_2=pti*sinh(eta);
    yi=0.5*log((yi_1+yi_2)/(yi_1-yi_2));

    // cout << yi << endl;

    x=sqrt(m*m+pti*pti)*exp(abs(yi)-yb)/mb; // lightcone variable

    cout << x << endl;

    Jacobian=sqrt(1-(m*m)/((m*m+pt*pt)*pow(cosh(y),2)));

    double R_eta; // Ridge yi에 대한 식을 eta로 바꾼 것
    
    if (x>=1) R_eta=0.;
    else R_eta=Jacobian*pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti); // Aridge을 포함하지 않은 eta에 대한 Ridge

    return (R_eta);
}