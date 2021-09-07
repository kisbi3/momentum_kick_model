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

//Aridge를 구하기 위해 적분할 함수
double integralAridge(double pti, double yi, int check){
    // std::cout<<pti<<std::setw(7)<<yi<<std::setw(9)<<phii<<std::endl;       
    double x = lightcone(pti, yi);
    double squareroot=sqrt(m*m+pti*pti);
    // std::cout<<pti<<std::setw(10)<<yi<<std::setw(10)<<phii<<std::setw(10)<<x<<std::endl;

    if(check == 1){
        if(x>=1.){
            return 0.;
        }
        else{
            return pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
        }
        
    }

    else{
        if(x>=1.){
            return 0.;
        }
        else{
            return pti*pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
        }
    }
    
}

//Ridge파트 적분
double RidgeDis(double Aridge, double ptf, double etaf, double phif, int check){
    double pti = sqrt(ptf*ptf-2*ptf*q*cos(phif)/cosh(etajet)+q*q/(cosh(etajet)*cosh(etajet)));
    // std::cout<<pti<<std::endl;
    double yi = rapidityintit(pti,etaf);
    double yf = rapidityintit(ptf,etaf);
    double E = sqrt(mb*mb+ptf*ptf)*cosh(yf);
    double Ei = sqrt(mb*mb+pti*pti)*cosh(yi);
    double x = lightcone(pti, yi);

    if (x>=1.){
        return 0.;
    }
    
    else{                                                                                                       // Ei/E 임을 명심하자.
        return (Aridge*integralAridge(pti, yi, check))*sqrt(1.-((mb*mb)/((mb*mb+ptf*ptf)*cosh(yf)*cosh(yf))))*Ei/E;
    }//
    
}

double Eii(double Aridge, double ptf, double etaf, double phif, int check){
    
    // double pti = ptf - q;
    // std::cout<<pti<<std::endl;

    if (check == 0){
        double pti = ptf-q;
        double yi = rapidityintit(pti,etaf);
        double Ei = sqrt(mb*mb+pti*pti)*cosh(yi);
        double x = lightcone(pti, yi);
        if(x>=1.){
            return 0.;
        }
        else{
            return Ei;
        }
    }
    
    else{       
        double pti = sqrt(ptf*ptf-2*ptf*q*cos(phif)/cosh(etajet)+q*q/(cosh(etajet)*cosh(etajet)));
        double yi = rapidityintit(pti,etaf);
        double yf = rapidityintit(ptf,etaf);
        double E = sqrt(mb*mb+ptf*ptf)*cosh(yf);
        double Ei = sqrt(mb*mb+pti*pti)*cosh(yi);
        double x = lightcone(pti, yi);
        if(x>=1.){
            return 0.;
        }
        else{
            return Ei;
        }
    }//
}

double Eff(double Aridge, double ptf, double etaf, double phif, int check){
    double pti = sqrt(ptf*ptf-2*ptf*q*cos(phif)/cosh(etajet)+q*q/(cosh(etajet)*cosh(etajet)));
    // double pti = ptf - q;
    // std::cout<<pti<<std::endl;
    double yi = rapidityintit(pti,etaf);
    double yf = rapidityintit(ptf,etaf);
    double E = sqrt(mb*mb+ptf*ptf)*cosh(yf);
    double Ei = sqrt(mb*mb+pti*pti)*cosh(yi);
    double x = lightcone(pti, yi);

    if (x>=1.){
        return 0.;
    }
    
    else{                                                                                                       // Ei/E 임을 명심하자.
        return E;
    }//
}

//Jetdistribution
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
//여기까지 Jet

int main()
{
    using std::cout;
    using std::endl;
    using std::setw;
    // std::string buffer;
    double dyi, dphii, sum, totalsum, phii, yi, dpti, pti, sum2, resultsum;
    int i, j, k, nyi, npti, nphii, check2;

    FILE *fpt;

    double Aridge=0.1420386444451514;

    fpt = fopen("Ei_E.csv", "w+");
    fprintf(fpt,"pt, Ei, Ef, Ei/Ef, Ef/Ei\n");

    double ptf, etaf, phif, detaf, dphif, ohno, etaf0, phif0, dist, Ei, Ef;
    int netaf, nphif;
    // cout<<"1"<<endl;

    netaf = 100;
    nphif = 100;

    detaf = double((1.4+1.4)/netaf);
    dphif = double((M_PI+M_PI)/nphif);

    sum = totalsum = 0.;
    etaf0 = -1.4;
    phif0 = -M_PI;

    check2 = 1;
    for(ptf=0.15;ptf<=4;ptf=ptf+0.01){
        phif = -M_PI;
        // std::cout<<ptf-q<<std::endl;
        for(j=1;j<=nphif;j+=1){
            etaf = -1.4;
            for (i=1;i<=netaf;i+=1){
                Ei += Eii(Aridge, ptf, etaf, phif, check2)*detaf*dphif;
                Ef += Eff(Aridge, ptf, etaf, phif, check2)*detaf*dphif;
                etaf += detaf;
            }
            
            // totalsum += sum*dphif;
            phif += dphif;

            
            // sum = 0.;

        }
        // cout<<ptf<<setw(20)<<sum<<endl;

        fprintf(fpt,"%f, %f, %f, %f, %f\n", ptf, Ei, Ef, Ei/Ef, Ef/Ei);

        Ei = Ef = 0.;
        // resultsum = resultsum + totalsum;
        // cout<<pt<<std::setw(20)<<totalsum<<endl;

    }



    fpt = fopen("Ei_E_phi.csv", "w+");
    fprintf(fpt,"pt, Ef, Ei_phi=0, Ei_phi=0.5, Ei_phi=1, Ei_phi=pi, Ei_phi=pi/2, Ei_phi=pi/4\n");

    double Ei_0, Ei_5, Ei_1, Ei_pi, Ei_pi2, Ei_pi4;

    detaf = double((1.4+1.4)/netaf);

    sum = totalsum = 0.;
    etaf0 = -1.4;
    phif0 = -M_PI;

    check2 = 1;
    for(ptf=0.15;ptf<=4;ptf=ptf+0.01){
        // std::cout<<ptf-q<<std::endl
        etaf = -1.4;
            for (i=1;i<=netaf;i+=1){
                Ei_0 += Eii(Aridge, ptf, etaf, 0., check2)*detaf;
                Ei_5 += Eii(Aridge, ptf, etaf, 0.5, check2)*detaf;
                Ei_1 += Eii(Aridge, ptf, etaf, 1., check2)*detaf;
                
                Ei_pi += Eii(Aridge, ptf, etaf, M_PI, check2)*detaf;
                Ei_pi2 += Eii(Aridge, ptf, etaf, M_PI/2, check2)*detaf;
                Ei_pi4 += Eii(Aridge, ptf, etaf, M_PI/4, check2)*detaf;
                Ef += Eff(Aridge, ptf, etaf, 0., check2)*detaf;


                etaf += detaf;
            }
            


        // cout<<ptf<<setw(20)<<sum<<endl;

        fprintf(fpt,"%f, %f, %f, %f, %f, %f, %f, %f\n", ptf, Ef, Ei_0, Ei_5, Ei_1, Ei_pi, Ei_pi2, Ei_pi4);

        Ei_0 = Ef = Ei_1 = Ei_5 = Ei_pi = Ei_pi2 = Ei_pi4 = 0.;
        // resultsum = resultsum + totalsum;
        // cout<<pt<<std::setw(20)<<totalsum<<endl;

    }



    fclose(fpt);



    return 0;
}
