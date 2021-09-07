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
double m = 0.13957018;  //m == mpi
double mb = m; //mb==mpi, GeV
// double md = m; //GeV -> 이는 md=mpi의 그래프를 그려보기 위함. 이때 Aridge를 따로 구하는게 맞는가?
double md = 1.;//GeV
double sqrSnn = 200.; //GeV
double mp = 0.938272046; //Proton mass, GeV
double Njet=0.75;
double fj=0.632;
double Tjet=0.55;
double fRNk = 4.;
double etajet = 0.;



double rapidityintit(double pt, double eta){
    double root=sqrt((pt*pt*cosh(eta)*cosh(eta))+mb*mb);
    double a = root+pt*sinh(eta);
    double b = root-pt*sinh(eta);
    return (log(a/b))/2;
}

double lightcone(double pti, double yi){
    // double yi = rapidityinit(pti);
    double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mpi
    double squareroot=sqrt(m*m+pti*pti);
    double yiabs = std::fabs(yi);
    // std::cout<<pti<<std::setw(8)<<yi<<std::setw(15)<<exp(yiabs-yb)<<std::setw(15)<<squareroot/mb<<std::endl;
    return (squareroot/m)*exp(yiabs-yb);
}

double ptii(double ptf, double phif){
    return sqrt(ptf*ptf-2*ptf*q*cos(phif)+q*q);
}

int main()
{
    using std::cout;
    using std::endl;
    using std::setw;
    // std::string buffer;
    double dyi, dphii, sum, totalsum, phii, yi, dpti, pti, sum2, resultsum;
    int i, j, k, nyi, npti, nphii, check2;

    FILE *fpt;

    fpt = fopen("pti,q=1.csv", "w+");
    fprintf(fpt,"ptf, pti, rapidity\n");

    double rapid;

    double ptf, etaf, phif, detaf, dphif, ohno, etaf0, phif0, dist, Ei, Ef;
    int netaf, nphif;
    // cout<<"1"<<endl;

    netaf = 100;
    nphif = 100;

    detaf = double((1.4+1.4)/netaf);
    dphif = double((1.+1.)/nphif);

    sum = totalsum = 0.;
    etaf0 = -1.4;
    phif0 = -1.;

    check2 = 1;
    for(ptf=0.15;ptf<=4;ptf=ptf+0.01){
        
        etaf = -1.4;
        // std::cout<<ptf-q<<std::endl;
        for(j=1;j<=netaf;j+=1){
            phif = -1.;
            for (i=1;i<=nphif;i+=1){
                pti += ptii(ptf, phif)*dphif*detaf;
                phif += dphif;
            }
            
            rapid += rapidityintit(ptf, etaf)*detaf;
            etaf += detaf;

        }

        fprintf(fpt,"%f, %f, %f\n", ptf, pti, rapid);

        pti = 0.;
        // resultsum = resultsum + totalsum;
        // cout<<pt<<std::setw(20)<<totalsum<<endl;

    }


    fpt = fopen("pti_3d,q=1.csv", "w+");
    fprintf(fpt,"ptf, phif, pti, pti/ptf\n");

    netaf = 100;
    nphif = 100;

    detaf = double((1.4+1.4)/netaf);
    dphif = double((M_PI+M_PI)/nphif);

    sum = totalsum = 0.;
    etaf0 = -1.4;
    phif0 = -1.;

    check2 = 1;
    for(ptf=0.15;ptf<=4;ptf=ptf+0.01){
        phif = -M_PI;
        // std::cout<<ptf-q<<std::endl;
        for(j=1;j<=nphif;j+=1){
            pti = ptii(ptf, phif);
        
            fprintf(fpt,"%f, %f, %f, %f\n", ptf, phif, pti, pti/ptf);
            phif += dphif;
            pti = 0.;

        }


        // cout<<ptf<<setw(20)<<sum<<endl;

        
        // resultsum = resultsum + totalsum;
        // cout<<pt<<std::setw(20)<<totalsum<<endl;

    }


    fpt = fopen("rapidity_3d.csv", "w+");
    fprintf(fpt,"ptf, etaf, rapidity\n");

    int n = 100;
    double dptf = double ((4-0.15)/n);
    detaf = double ((1.4+1.4)/n);
    dphif = double ((1.+1.)/n);
    etaf = -1.4;
    phif = -1.;

    for(j=1;j<=n+1;j+=1){
        ptf = 0.15;
        for (i=1;i<=n;i+=1){  
            fprintf(fpt, "%f, %f, %f\n", etaf, ptf, rapidityintit(ptf, etaf));
            ptf += dptf;
        }
        etaf += detaf;
    }




    fclose(fpt);



    return 0;
}
