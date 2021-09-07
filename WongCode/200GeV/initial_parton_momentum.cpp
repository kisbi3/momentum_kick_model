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
    // double pti = sqrt(ptf*ptf-2*ptf*q*cos(phif)/cosh(etajet)+q*q/(cosh(etajet)*cosh(etajet)));
    // std::cout<<pti<<std::endl;
    // double yi = rapidityintit(ptf,etaf);
    double yf = rapidityintit(ptf,etaf);
    // double E = sqrt(mb*mb+ptf*ptf)*cosh(yf);
    // double Ei = sqrt(mb*mb+pti*pti)*cosh(yi);
    double x = lightcone(ptf, yf);

    if (x>=1.){
        return 0.;
    }
    
    else{
        return (Aridge*integralAridge(ptf, yf, check))*sqrt(1.-((mb*mb)/((mb*mb+ptf*ptf)*cosh(yf)*cosh(yf))));
    }//
    
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
    // fpt = fopen("md=mpi_Ridge.csv", "w+");

    
    //Aridge를 구하기 위한 적분

    //pti -> 0~5, yi -> -6~+6, phii -> 0~2pi

    nyi = 10000;
    npti = 10000;
    nphii = 1000;

    dyi = double((0.0+10.0)/nyi);
    dpti = double((0.+20.0)/npti);
    dphii = double (M_PI+M_PI)/nphii;

    sum = totalsum = 0.;
    pti = 0.0;  //적분을 pt=0부터 하는것이 옳은가? 원통좌표계에서의 적분인데?
    yi = 0.;    //0~4 적분한 후 x2할 것.
    double yi0 = 0.;
    
    double check = 0.;
    double intaridge = 0.;

    // cout<<dyi<<endl;

    check2 = 0;
    for(k=1;k<=npti+1;k+=1){
        for (i=1;i<=nyi;i+=1){

            yi += dyi;
            check = intaridge;
            intaridge = integralAridge(pti, yi, check2);

            // cout<<yi<<endl;

            sum += intaridge;

            if(i == nyi){
                // cout<<"??"<<endl;
                sum += (integralAridge(pti, yi0, check2)-intaridge)/2.;
                sum *= dyi;
                break;
            }
            else if(i != 1 && intaridge == 0.){
                // cout<<"!!!"<<endl;
                // sum -= check;
                sum += (integralAridge(pti, yi0, check2)-check)/2.;
                sum *= dyi;
                break;       
            }

        }
        yi = 0.;
        sum *= 2.;
        // cout<<sum<<endl;
        // sum -= check;
        double checksum;
        totalsum += sum;
        // cout<<pti<<endl;
        if (k==1){
            checksum = totalsum;
        }

        else if (k==npti+1 || lightcone(pti+dpti,0.)>=1.){
            
            totalsum -= (sum+checksum)/2.;
            // cout<<'2'<<endl;
            // cout<<pti<<setw(15)<<sum<<setw(15)<<totalsum<<endl;
            totalsum *= dpti;
            break;
            
        }
        // cout<<pti<<setw(15)<<sum<<setw(15)<<totalsum<<endl;
        
        pti += dpti;
        sum = 0.;

        // cout<<pt<<std::setw(20)<<totalsum<<endl;
    }



    totalsum *= 2*M_PI;

    // cout<<lightcone(pti+dpti,0.0)<<setw(15)<<integralAridge(pti+dpti,0.)<<endl;

    // fclose(fpt);
    cout.precision(10);
    double Aridge = 1/totalsum;
    // cout<<totalsum<<setw(15)<<Aridge<<endl;
    
    // double Aridge = 1/resultsum;
    cout<<totalsum<<setw(15)<<Aridge<<endl;
    
    //Aridge=0.1420386444451514

    // Aridge=0.1420386444451514;
    //여기까지 Aridge를 구하기 위한 적분




    //여기부터 Ridge항 적분 pt : 0.15~4, etaf : -1.4~+1.4, phif : -1~1

    fpt = fopen("Initial_parton.csv", "w+");
    fprintf(fpt,"pt, RidgeDis\n");

    double ptf, etaf, phif, detaf, dphif, ohno, etaf0, phif0, dist;
    int netaf, nphif;
    cout<<"1"<<endl;

    netaf = 100;
    nphif = 100;

    detaf = double((1.4+1.4)/netaf);
    dphif = double((1.0+1.0)/nphif);

    sum = totalsum = 0.;
    etaf0 = -1.4;
    phif0 = -1.0;

    check2 = 1;
    for(ptf=0.15;ptf<=4;ptf=ptf+0.01){
        phif = -1.;
        // std::cout<<ptf-q<<std::endl;
        for(j=1;j<=nphif;j+=1){
            etaf = -1.4;
            for (i=1;i<=netaf;i+=1){
                
                sum += RidgeDis(Aridge, ptf, etaf, phif, check2)*detaf*dphif;
                etaf += detaf;
            }
            
            // totalsum += sum*dphif;
            phif += dphif;

            
            // sum = 0.;

        }
        // cout<<ptf<<setw(20)<<sum<<endl;

        sum *= 8/3;

        fprintf(fpt,"%f, %f\n", ptf, sum);

        sum = 0.;
        // resultsum = resultsum + totalsum;
        // cout<<pt<<std::setw(20)<<totalsum<<endl;

    }


    fclose(fpt);



    return 0;
}
