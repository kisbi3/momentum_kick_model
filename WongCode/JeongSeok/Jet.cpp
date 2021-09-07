#include <iostream>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#define PI 3.1415926535897932385
#define mb 0.938272046 // proton mass
#define m 0.13957018 // pion mass

double Njetdis(double eta, double phi, double pt);
double Ieta(double pt, double phi, double eta0, double etan, int n);
double Iphi(double pt, double phi0, double phin, double eta0, double etan, int n);

int main()
{
    FILE *fpt;

    fpt = fopen("Jet.csv", "w+");
    fprintf(fpt,"pt, Njetdis\n");
    
    int n=1000;
    double pt;
    double eta0, etan, phi0, phin, fJ, integral, Njdis;
    eta0=-1.4;
    etan=1.4;
    phi0=-1.0;
    phin=1.0;
    // fJ=0.632;

    for(pt=0.15;pt<=4;pt=pt+0.01)
    {
        Njdis=Iphi(pt,phi0,phin,eta0,etan,n);
        fprintf(fpt,"%f, %f\n", pt, Njdis);
    }
    fclose(fpt);
}

double Njetdis(double eta, double phi, double pt)
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

double Ieta(double pt, double phi, double eta0, double etan, int n)
{
    double eta=0, sum_mid=0;
    double deta=(etan-eta0)/n;
    int i;

    for(i=1;i<n;i++)
    {
        eta=eta0+i*deta;
        sum_mid+=Njetdis(eta,phi,pt);
    }
    return ((Njetdis(eta0,phi,pt)+Njetdis(etan,phi,pt)+2.*sum_mid)*deta/2.);
}

double Iphi(double pt, double phi0, double phin, double eta0, double etan, int n)
{
    int i;
    double phi,sum_mid=0.,dphi;

    dphi=(phin-phi0)/n;

    for(i=1;i<n;i++)
    {
        phi = phi0+i*dphi;
        sum_mid += Ieta(pt,phi,eta0,etan,n);
    }
    return((Ieta(pt,phi0,eta0,etan,n)+Ieta(pt,phin,eta0,etan,n)+2.*sum_mid)*dphi/2.);
}