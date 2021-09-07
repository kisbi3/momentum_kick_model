#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstring>

//JF = Jet Fragments, MP = Medium Parton

double pt_trig = 5.; //GeV/c
double T_JF = 0.266+0.084*pt_trig;  //GeV (43)
double N_JF = 0.15+0.1*pt_trig;
double ma = 1.1;    //GeV
double sigmaphizero=0.55;   //논문에 언급이 없는데?
double m=0.1396; // m=mpi(파이온 질량)
// double Njet=0.75;
// double Tjet=0.55;

//sigmaphi square
double sigmaphisq(double pt){
    return (sigmaphizero*sigmaphizero*ma*ma)/(ma*ma+pt*pt);
}

double integralNjet(double pt, double eta, double phi, double constant){
    double sigmaphi;
    sigmaphi = sigmaphisq(pt);
    return (constant/sigmaphi)*exp(((m-sqrt(m*m+pt*pt))/T_JF)-((phi*phi+eta*eta)/(2*sigmaphi)));
}

int main()
{
    using std::cout;
    using std::endl;
    using std::setw;
    // std::string buffer;
    double pt, eta, phi, deta, dphi, constant, sum, totalsum;
    int i, j, nphi, neta;

    FILE *fpt;
    fpt = fopen("7TeV_Jet.csv", "w+");
    fprintf(fpt,"pt, Njetdis\n");

    neta = 100;
    nphi = 100;

    dphi = double((1.5+1.5)/nphi);
    deta = double((4.8+4.8)/neta);


    //dphi = dx
    //deta = dy
    eta = -4.8;
    phi = -1.5;
    constant = N_JF/(T_JF*(m+T_JF)*2*M_PI);
    double intnjet = 0.;

    cout<<N_JF<<setw(15)<<T_JF<<endl;

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
            phi = -1.5;
            eta += deta;
        }
        eta = -4.8;
        // cout<<pt<<std::setw(20)<<totalsum<<endl;
        
        fprintf(fpt,"%f, %f\n", pt, totalsum);
        totalsum = 0.;
    }

    fclose(fpt);



    return 0;
}
