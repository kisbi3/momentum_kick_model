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
double m = 0.13957018;  //m == mpi
double mb = m; //mb==mpi, GeV
double sqrSnn = 200; //GeV
double mp = 0.938; //Proton mass, GeV


//Aridge를 구하기 위해 적분할 함수
double integralAridge(double pti, double yi){
    return pti*pti*pti+yi*yi;
}

int main()
{
    using std::cout;
    using std::endl;
    using std::setw;
    // std::string buffer;
    double dyi, dphii, sum, totalsum, phii, yi, dpti, pti;
    int i, j, k, nyi, npti;
    
    //Aridge를 구하기 위한 적분

    //pti -> 0~5, yi -> -6~+6, phii -> 0~2pi

    nyi = 10;
    npti = 10;

    dyi = double((0.0+1.0)/nyi);
    dpti = double((0.0+1.0)/npti);

    sum = totalsum = 0.;
    pti = 0.0;
    yi = 0.0;
    
    for(k=1;k<=npti+1;k+=1){
        for (i=1;i<=nyi+1;i+=1){
            cout<<pti<<setw(10)<<yi<<setw(10)<<integralAridge(pti, yi)<<endl;
            if (i==1){
                sum += integralAridge(pti, yi)/2.;
            }
            else if (i==nyi+1){
                sum += integralAridge(pti, yi)/2.;
                sum = sum*dyi;
            }
            else {
                sum += integralAridge(pti, yi);
            }
            
            yi += dyi;
        }
        if (k==1){
            totalsum += sum/2.;
        }
        else if (k==npti+1){
            totalsum += sum/2.;
            totalsum = totalsum*dpti;
        }
        else{
            totalsum += sum;
        }

        yi = 0.0;
        pti += dpti;
        sum = 0.;

        // resultsum = resultsum + totalsum;
        // cout<<pt<<std::setw(20)<<totalsum<<endl;
        
        // fprintf(fpt,"%f, %f\n", pt, totalsum);
        // totalsum = 0.;
    }

    // fclose(fpt);
    double Aridge = 1/totalsum;
    cout<<totalsum<<setw(10)<<Aridge<<endl;


    return 0;
}
