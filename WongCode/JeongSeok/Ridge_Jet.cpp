#include <iostream>
#include <stdio.h>
#include <math.h>
#define PI 3.1415926535897932385
#define mb 0.938272046 // proton mass
#define m 0.13957018 // pion mass
// using namespace std;

double Jet(double eta, double phi, double pt);
double IJeta(double pt, double phi, double eta0, double etan, int n);
double IJphi(double pt, double phi0, double phin, double eta0, double etan, int n);
double Ridge(double pt, double eta);
double IReta(double pt, double eta0, double etan, int n);

int main()
{
    FILE *fpt;

    fpt = fopen("Ridge_Jet_asdf.csv", "w+");
    fprintf(fpt,"pt, N_Rdige_dis\n");

    int n=1000;
    double pt, eta0, etan, phi0, phin, Ridge_Jet;
    eta0=-1.4;
    etan=1.4;
    phi0=-1.0;
    phin=1.0;

    for(pt=0.15;pt<=4;pt=pt+0.01)
    {
        // Ridge_Jet=(8./3.)*IReta(pt,eta0,etan,n)+0.632*IJphi(pt,phi0,phin,eta0,etan,n);
        Ridge_Jet=IReta(pt,eta0,etan,n);
        fprintf(fpt,"%f, %f\n",pt,Ridge_Jet);
    }
    fclose(fpt);
}

double Jet(double eta, double phi, double pt)
{
    double Njet, Tjet;
    Njet=0.75;
    Tjet=0.55;

    double sigmaphi0, ma, sigmaphi;
    sigmaphi0=0.5;
    ma=1.1;
    sigmaphi=(sigmaphi0*ma)/sqrt(ma*ma+pt*pt);

    return ((Njet*exp((m-sqrt(m*m+pt*pt))/Tjet)/(Tjet*(m+Tjet)))*(exp(-(phi*phi+eta*eta)/(2.*sigmaphi*sigmaphi))/(2.*PI*sigmaphi*sigmaphi)));
}

double Ridge(double pt, double eta)
{
    double T, md, a, q, sqrt_sNN, yb, Aridge; // constant
    T=0.70;
    md=1.0;
    a=0.5;
    q=2.0;
    sqrt_sNN=7000.;
    yb=acosh(sqrt_sNN/(2.*mb));
    Aridge = 0.0399224765;
    //내 값으로 일단 변경
    
    double pti, y_1, y_2, y, yi_1, yi_2, yi, x, Jacobian, E, Ei; // parameter
    pti=pt-q;
    y_1=sqrt(pow(pt*cosh(eta),2)+m*m);
    y_2=pt*sinh(eta);
    y=0.5*log((y_1+y_2)/(y_1-y_2));
    
    yi_1=sqrt(pow(pti*cosh(eta),2)+m*m);
    yi_2=pti*sinh(eta);
    yi=0.5*log((yi_1+yi_2)/(yi_1-yi_2));

    x=sqrt(m*m+pti*pti)*exp(abs(yi)-yb)/mb; // lightcone variable
    Jacobian=sqrt(1-(m*m)/((m*m+pt*pt)*pow(cosh(y),2)));
    E=sqrt(m*m+pt*pt)*cosh(y);
    Ei=sqrt(m*m+pti*pti)*cosh(yi);

    double R_eta; // Ridge yi에 대한 식을 eta로 바꾼 것
    if (x>1) R_eta=0.;
    else R_eta=Aridge*Jacobian*(E/Ei)*pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);; // Aridge을 포함하지 않은 eta에 대한 Ridge

    return (R_eta);
}


double IReta(double pt, double eta0, double etan, int n)
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
    Integral=2.*PI*(Ridge(pt,eta0)+Ridge(pt,etan)+2.*sum_mid)*deta/2.;

    return (Integral);
}

double IJeta(double pt, double phi, double eta0, double etan, int n)
{
    double eta=0, sum_mid=0;
    double deta=(etan-eta0)/n;
    int i;

    for(i=1;i<n;i++)
    {
        eta=eta0+i*deta;
        sum_mid+=Jet(eta,phi,pt);
    }
    return ((Jet(eta0,phi,pt)+Jet(etan,phi,pt)+2.*sum_mid)*deta/2.);
}

double IJphi(double pt, double phi0, double phin, double eta0, double etan, int n)
{
    int i;
    double phi,sum_mid=0.,dphi;

    dphi=(phin-phi0)/n;

    for(i=1;i<=n-1;i++)
    {
        phi = phi0+i*dphi;
        sum_mid += IJeta(pt,phi,eta0,etan,n);
    }
    return((IJeta(pt,phi0,eta0,etan,n)+IJeta(pt,phin,eta0,etan,n)+2.*sum_mid)*dphi/2.);
}