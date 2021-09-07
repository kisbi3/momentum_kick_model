#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstring>

//sigmaphi square
double sigmaphisq(double pt){
    double sigmaphizero, ma;
    sigmaphizero = 0.5;
    ma = 1.1;
    return (sigmaphizero*sigmaphizero*ma*ma)/(ma*ma+pt*pt);
    
}

double integralNjet(double pt, double eta, double phi, double constant){
    double sigmaphi, Tjet, m;
    sigmaphi = sigmaphisq(pt);
    Tjet=0.55;
    m=0.1396;
    return (constant/sigmaphi)*exp(((m-sqrt(m*m+pt*pt))/Tjet)-((phi*phi+eta*eta)/(2*sigmaphi)));
}

int main()
{
    using std::cout;
    using std::endl;
    using std::setw;
    // std::string buffer;
    double Njet, fj, sigma0, ma, pt, eta, phi, deta, dphi, m, constant, Tjet, sum, totalsum;
    int i, j, nphi, neta;

    FILE *fpt;
    fpt = fopen("Jet.csv", "w+");
    fprintf(fpt,"pt, Njetdis\n");

    neta = 1000;
    nphi = 1000;

    dphi = double((1.0+1.0)/nphi);
    deta = double((1.4+1.4)/neta);


    Njet=0.75;
    fj=0.632;
    sigma0=0.5;
    ma=1.1;
    m=0.1396; // m=mpi(파이온 질량)
    Tjet=0.55;

    //dphi = dx
    //deta = dy
    eta = -1.4;
    phi = -1.;
    constant = Njet/(Tjet*(m+Tjet)*2*M_PI);
    double intnjet = 0.;

    for(pt=0.15;pt<=4;pt=pt+0.01){
        for (i=1;i<=neta+1;i+=1){
            for (j=1;j<=nphi+1;j+=1){
                intnjet = integralNjet(pt, eta, phi, constant);
                if (j==1){
                    sum += intnjet/2.;
                }
                else if (j==nphi+1){
                    sum += intnjet/2.;
                    sum = sum*dphi;
                }
                else {
                    sum += intnjet;
                }
                // std::cout<<pt<<std::setw(15)<<eta<<std::setw(15)<<phi<<std::setw(15)<<intnjet<<std::endl;
                phi += dphi;
            }
            

            if (i==1){
                    totalsum += sum/2.;
                }
            else if (i==neta+1){
                totalsum += sum/2.;
                totalsum = totalsum*deta;
            }
            else{
                totalsum += sum;
            }
            sum = 0.;
            phi = -1.;
            eta += deta;
        }
        eta = -1.4;
        // cout<<pt<<std::setw(20)<<totalsum<<endl;
        
        fprintf(fpt,"%f, %f\n", pt, totalsum);
        totalsum = 0.;
    }

    fclose(fpt);



    return 0;
}
