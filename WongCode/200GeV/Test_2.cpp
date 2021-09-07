#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstring>

double a = 0.5; //fall off parameter
double T = 0.5; //Temperature, Gev
double q = 1.;  //GeV
double fRN = 4.;
double md = 1.; //GeV
double mb = 0.13957018; //mb==mpi, GeV
double sqrSnn = 200; //GeV
double mp = 0.938; //Proton mass, GeV


// sigmaphi square
// double sigmaphisq(double pt){
//     double sigmaphizero, ma;
//     sigmaphizero = 0.5;
//     ma = 1.1;
//     return (sigmaphizero*sigmaphizero*ma*ma)/(ma*ma+pt*pt);
    
// }
// double rapidityintit(double pti){
//     double root=sqrt((pti*pti*cosh(eta)*cosh(eta))+mb*mb);
//     double a = root+pti*sinh(eta);
//     double b = root-pti*sinh(eta);
//     return (log(a/b))/2;
// }

double lightcone(double pti, double yi){
    // double yi = rapidityinit(pti);
    double yb = acosh(sqrSnn/(2.*mb));    //mN=mbeam, mb = mpi
    double squareroot=sqrt(mb*mb+pti*pti);
    double yiabs = std::fabs(yi);
    // std::cout<<pti<<std::setw(8)<<yi<<std::setw(15)<<exp(yiabs-yb)<<std::setw(15)<<squareroot/mb<<std::endl;
    return (squareroot/mb)*exp(yiabs-yb);
}

double integralAridge(double pti, double yi, double phii){
    // std::cout<<pti<<std::setw(7)<<yi<<std::setw(9)<<phii<<std::endl;       
    double x = lightcone(pti, yi);
    double squareroot=sqrt(mb*mb+pti*pti);
    // std::cout<<pti<<std::setw(10)<<yi<<std::setw(10)<<phii<<std::setw(10)<<x<<std::endl;
    if (x>=1.){
        return 0.;
    }
    return pti*pow((1-x),a)*((exp(-squareroot/T))/(sqrt(md*md+(pti-1)*(pti-1))));
}

int main()
{
    using std::cout;
    using std::endl;
    using std::setw;
    // std::string buffer;
    double pt, dyi, dphii, sum, totalsum, phii, yi, dpti, pti, resultsum;
    int i, j, k, nphii, nyi, npti;

    

    // FILE *fpt;
    // fpt = fopen("Me.csv", "w+");
    // fprintf(fpt,"pt, RidgeDis\n");

    nyi = 100;
    nphii = 100;
    npti = 100;

    dyi = double((6.0+6.0)/nyi);
    dpti = double((0.0+5.0)/npti);
    dphii = double((0+2*M_PI)/nphii);

    sum = resultsum = totalsum = 0.;
    pti = 0.;
    yi = -6.;
    phii = 0.;

    for(k=1;k<=nphii;k+=1){
        for (i=1;i<=nyi;i+=1){
            for (j=1;j<=npti;j+=1){
                // cout<<integralAridge(pt, yi, phii)<<setw(10);
                if (j==1){
                    sum += integralAridge(pti, yi, phii)/2.;
                }
                else if (j==npti){
                    sum += integralAridge(pti, yi, phii)/2.;
                    sum = sum*dpti;
                }
                else {
                    sum += integralAridge(pti, yi, phii);
                }
                
                // cout<<pti<<std::setw(7)<<yi<<std::setw(9)<<phii<<endl;       
                // cout<<mb*mb+pti*pti<<setw(10)<<sqrt(mb*mb+pti*pti)<<endl;
                pti += dpti;
            }
           



            if (i==1){
                    totalsum += sum/2.;
                }
            else if (i==nyi){
                totalsum += sum/2.;
                totalsum = totalsum*dyi;
            }
            else{
                totalsum += sum;
            }
            sum = 0.;
            pti = 0.;
            yi += dyi;
        }

        if (k==1){
            resultsum += totalsum/2.;
        }
        else if (k==nphii){
            resultsum += totalsum/2.;
            resultsum = resultsum*dphii;
        }
        else{
            resultsum += totalsum;
        }

        yi = -6.;
        phii += dphii;
        totalsum = 0.;

        // resultsum = resultsum + totalsum;
        // cout<<pt<<std::setw(20)<<totalsum<<endl;
        
        // fprintf(fpt,"%f, %f\n", pt, totalsum);
        // totalsum = 0.;
    }

    // fclose(fpt);
    double Aridge = 1/resultsum;
    cout<<resultsum<<setw(10)<<Aridge<<endl;

    return 0;
}
