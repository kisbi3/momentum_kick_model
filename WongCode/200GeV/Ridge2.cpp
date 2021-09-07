#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <fstream>

double a = 0.5; //fall off parameter
double T = 0.5; //Temperature, Gev
// double T = 1.3;
// double T = 1.;
double q = 1.;  //GeV
// double q = 1.2;
double m = 0.13957018;  //m == mpi
// double m = 0.14;  //m == mpi
double mb = m; //mb==mpi, GeV
// double md = m; //GeV -> 이는 md=mpi의 그래프를 그려보기 위함. 이때 Aridge를 따로 구하는게 맞는가?
double md = .77;//GeV
double sqrSnn = 200.; //GeV
// double sqrSnn = 13000.;
// double mp = 0.938272046; //Proton mass, GeV
double mp = 0.93890595;
double Njet=0.75;
double fj=0.632;
double Tjet=0.55;
double fRNk = 4.;
// double fRNk = 0.1;
// double fRNk = 1.;
// double fRNk1_2 = 1.0;
// double fRNk2_3 = 1.0;
// double fRNk3_4 = 1.0;
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
    double sigmay = 6.6;

    if(check == 1){
        // if(x>=1.){
        //     return 0.;
        // }
        // else{
        //     return pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
        // }
        return exp(-yi*yi/(2*sigmay*sigmay))*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
        
    }

    //Aridge
    else{
        // if(x>=1.){
        //     return 0.;
        // }
        // else{
        //     // return pti*pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
        
        // }
        return pti*exp(-yi*yi/(2*sigmay*sigmay))*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
    }
    
}

//E/Ei 확인
double E_Ei(double Aridge, double ptf, double etaf, double phif){
    double pti = sqrt(ptf*ptf-2*ptf*q*cos(phif)/cosh(etajet)+q*q/(cosh(etajet)*cosh(etajet)));
    double yi = rapidityintit(pti,etaf);
    double yf = rapidityintit(ptf,etaf);
    // double E = sqrt(ptf*ptf*cosh(etaf)*cosh(etaf)+m*m);
    // double Ei = sqrt(pti*pti*cosh(etaf)*cosh(etaf)+m*m);
    double E = sqrt(mb*mb+ptf*ptf)*cosh(yf);
    // double Ei = sqrt(mb*mb+pti*pti)*cosh(yi);
    double Ei = sqrt(pti*pti+ptf*ptf*sinh(etaf)*sinh(etaf)+m*m);

    return E/Ei;
}

//Ridge파트 적분
double RidgeDis(double Aridge, double ptf, double etaf, double phif, int check){
    double pti = sqrt(ptf*ptf-2*ptf*q*cos(phif)/cosh(etajet)+q*q/(cosh(etajet)*cosh(etajet)));
    // std::cout<<pti<<std::endl;
    // double yi = rapidityintit(pti,etaf);
    // double yf = rapidityintit(ptf,etaf);
    double E = sqrt(ptf*ptf*cosh(etaf)*cosh(etaf)+m*m);
    // double Ei = sqrt(pti*pti*cosh(etaf)*cosh(etaf)+m*m);
    double Ei = sqrt(pti*pti+ptf*ptf*sinh(etaf)*sinh(etaf)+m*m);
    // double E = sqrt(ptf*ptf+m*m);
    // double Ei = sqrt(pti*pti+m*m);
    // double E = sqrt(mb*mb+ptf*ptf)*cosh(yf);
    // double Ei = sqrt(mb*mb+pti*pti)*cosh(yi);

    double yi = log((Ei+ptf*sinh(etaf))/(Ei-ptf*sinh(etaf)))/2;
    double yf = log((E+ptf*sinh(etaf))/(E-ptf*sinh(etaf)))/2;

    double x = lightcone(pti, yi);

    double sigmay = 6.6;
    Aridge = 0.08;

    if (x>=1.){
        return 0.;
    }
    
    else{                                                                                                       // E/Ei 임을 명심하자.
        return (Aridge*exp(-yi*yi/(2*sigmay*sigmay))*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti))*sqrt(1.-((mb*mb)/((mb*mb+ptf*ptf)*cosh(yf)*cosh(yf))))*(E/Ei);
        // return (Aridge*integralAridge(pti, yi, check))*sqrt(1.-((mb*mb)/((mb*mb+ptf*ptf)*cosh(yf)*cosh(yf))))*(E/Ei);
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
    // fpt = fopen("md=mpi_Ridge.csv", "w+");

    
    //Aridge를 구하기 위한 적분

    //pti -> 0~5, yi -> -6~+6, phii -> 0~2pi

    nyi = 1000;
    npti = 1000;
    nphii = 1000;

    dyi = double((0.0+10.)/nyi);
    dpti = double((0.+10.)/npti);
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

    // for(k=1;k<=npti+1;k+=1){
    //     yi = 0.;
    //     for (i=1;i<=nyi;i+=1){
    //         totalsum += 2*integralAridge(pti, yi, check2)*dpti*dyi;
    //         yi += dyi;
    //         // cout<<pti<<setw(15)<<yi<<setw(15)<<integralAridge(pti, yi, check2)*dpti*dyi<<setw(20)<<totalsum<<endl;
    //     }
    //     // totalsum *= 2.;
    //     pti += dpti;
    // }


    totalsum *= 2*M_PI;

    // cout<<lightcone(pti+dpti,0.0)<<setw(15)<<integralAridge(pti+dpti,0.)<<endl;

    // fclose(fpt);
    cout.precision(10);
    double Aridge = 1/totalsum;
    // cout<<totalsum<<setw(15)<<Aridge<<endl;
    
    // double Aridge = 1/resultsum;
    cout<<totalsum<<setw(20)<<Aridge<<endl;
    
    //Aridge=0.1420386444451514

    // Aridge=0.1420386444451514;
    //여기까지 Aridge를 구하기 위한 적분


    

    //여기부터 Ridge항 적분 pt : 0.15~4, etaf : -1.4~+1.4, phif : -1~1

    fpt = fopen("Ridge2_phi.csv", "w+");
    fprintf(fpt,"pt, RidgeDis\n");
    double ptf, etaf, phif, detaf, dphif, ohno, etaf0, phif0, dist;
    int netaf, nphif;
    // cout<<"1"<<endl;

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
                
                // etaf += detaf;
                // double yi = rapidityintit(sqrt(ptf*ptf-2*ptf*q*cos(phif)+q*q),etaf+detaf);
                // check = dist;
                // dist = RidgeDis(Aridge, ptf, etaf, phif, check2);

                // // cout<<etaf<<endl;
                
                // sum += dist;
                sum += RidgeDis(Aridge, ptf, etaf, phif, check2)*detaf*dphif;
                etaf += detaf;
            }
            
            // totalsum += sum*dphif;
            phif += dphif;

            
            // sum = 0.;

        }
        // cout<<ptf<<setw(20)<<sum<<endl;

        sum = (sum*2*fRNk)/3;

        fprintf(fpt,"%f, %f\n", ptf, sum);

        sum = 0.;
        // resultsum = resultsum + totalsum;
        // cout<<pt<<std::setw(20)<<totalsum<<endl;

    }

    int n;
    /*

    // fpt = fopen("dNdy.csv", "w+");
    // fprintf(fpt,"y, dFdy\n");
    // n = 888;
    // dpti = double ((8.8-0.0)/n);
    // sum = 0.;


    // pti = 0.;
    // check2 = 0;

    // for(yi=0.;yi<=5.4;yi+=0.1){
    //     pti = 0.;
    //     for(j=1;j<=n;j+=1){
            
    //         sum += Aridge*integralAridge(pti, yi, check2)*dpti;
    //         pti += dpti;
    //     }
    //     fprintf(fpt,"%f, %f\n", yi, sum);
    //     sum = 0.;
    // }



    // fpt = fopen("Correlation.csv", "w+");
    // fprintf(fpt,"phi, dN/dphi\n");
    using std::ofstream;

    ofstream fout("/home/jaesung/Desktop/Dropbox/Code/WongCode/Correlation/phiCorrelation_pt1-2.csv");
    fout<<"phi,dN/dphi\n";
    double dptf;

    n = 100;
    dptf = double ((2.0-1.0)/n);
    detaf = double ((4.-2.)/n);
    sum = 0.;
    ptf = 0.;
    phif = -1.18;

    check2 = 1;

    for(phif=-1.18;phif<=1.18;phif+=0.005){
        etaf = 2.;
        // cout<<phif<<endl;
        for(i=1;i<=n;i++){
            ptf = 1.;
            // cout<<etaf<<endl;
            for(j=1;j<=n;j+=1){
                sum += ptf*RidgeDis(Aridge, ptf, etaf, phif, check2)*dptf*detaf;
                // cout<<phif<<setw(15)<<etaf<<setw(15)<<ptf<<setw(15)<<sum<<endl;
                ptf += dptf;
            }
            etaf += detaf;
        }
        sum *= 2*fRNk1_2/3;
        // fprintf(fpt,"%f, %f\n", phif, sum);

        fout<<phif<<","<<sum<<endl;
        sum = 0.;
    }
    // fout.close();



    ofstream fout23 ("/home/jaesung/Desktop/Dropbox/Code/WongCode/Correlation/phiCorrelation_pt2-3.csv");
    fout23<<"phi,dN/dphi\n";
    // double dptf;

    n = 100;
    dptf = double ((3.0-2.0)/n);
    detaf = double ((4.-2.)/n);
    sum = 0.;
    ptf = 0.;
    phif = -1.18;

    check2 = 1;

    for(phif=-1.18;phif<=1.18;phif+=0.005){
        etaf = 2.;
        // cout<<phif<<endl;
        for(i=1;i<=n;i++){
            ptf = 2.;
            // cout<<etaf<<endl;
            for(j=1;j<=n;j+=1){
                sum += ptf*RidgeDis(Aridge, ptf, etaf, phif, check2)*dptf*detaf;
                // cout<<phif<<setw(15)<<etaf<<setw(15)<<ptf<<setw(15)<<sum<<endl;
                ptf += dptf;
            }
            etaf += detaf;
        }
        sum *= 2*fRNk2_3/3;
        // fprintf(fpt,"%f, %f\n", phif, sum);

        fout23<<phif<<","<<sum<<endl;
        sum = 0.;
    }
    fout23.close();

    std::ofstream fout34("/home/jaesung/Desktop/Dropbox/Code/WongCode/Correlation/phiCorrelation_pt3-4.csv");
    fout34<<"phi,dN/dphi\n";
    // double dptf;

    n = 100;
    dptf = double ((4.0-3.0)/n);
    detaf = double ((4.-2.)/n);
    sum = 0.;
    ptf = 0.;
    phif = -1.18;

    check2 = 1;

    for(phif=-1.18;phif<=1.18;phif+=0.005){
        etaf = 2.;
        // cout<<phif<<endl;
        for(i=1;i<=n;i++){
            ptf = 3.;
            // cout<<etaf<<endl;
            for(j=1;j<=n;j+=1){
                sum += ptf*RidgeDis(Aridge, ptf, etaf, phif, check2)*dptf*detaf;
                // cout<<phif<<setw(15)<<etaf<<setw(15)<<ptf<<setw(15)<<sum<<endl;
                ptf += dptf;
            }
            etaf += detaf;
        }
        sum *= 2*fRNk3_4/3;
        // fprintf(fpt,"%f, %f\n", phif, sum);

        fout34<<phif<<","<<sum<<endl;
        sum = 0.;
    }
    fout34.close();


    ofstream fout14("/home/jaesung/Desktop/Dropbox/Code/WongCode/Correlation/phiCorrelation_pt1-4.csv");
    fout14<<"phi,dN/dphi\n";
    // double dptf;

    n = 100;
    dptf = double ((4.0-1.0)/n);
    detaf = double ((4.-2.)/n);
    sum = 0.;
    ptf = 0.;
    phif = -1.18;

    check2 = 1;

    for(phif=-1.18;phif<=1.18;phif+=0.005){
        etaf = 2.;
        // cout<<phif<<endl;
        for(i=1;i<=n;i++){
            ptf = 1.;
            // cout<<etaf<<endl;
            for(j=1;j<=n;j+=1){
                sum += ptf*RidgeDis(Aridge, ptf, etaf, phif, check2)*dptf*detaf;
                // cout<<phif<<setw(15)<<etaf<<setw(15)<<ptf<<setw(15)<<sum<<endl;
                ptf += dptf;
            }
            etaf += detaf;
        }
        sum *= 2*fRNk3_4/3;
        // fprintf(fpt,"%f, %f\n", phif, sum);

        fout14<<phif<<","<<sum<<endl;
        sum = 0.;
    }
    fout14.close();


*/




/*

    //E/Ei 3차원 그래프 그리기

    // double dptf, plots;

    std::ofstream fout2("E_Ei/E_Ei_3d.csv");
    fout2<<"eta,phi,E/Ei\n";
    // double result[10][1000] = {0.};
    

    // double detaf, dphif, etaf, phif, ptf;

    n = 100;
    dptf = double ((2.-1.)/n);
    detaf = double ((5.+5.)/n);
    dphif = double ((1.+1.)/n);
    etaf = -5.;
    phif = -1.;
    double ptf0 = 1.;
    sum = 0.;


    for(j=1;j<=n+1;j+=1){
        phif = -1.;
        for(k=1;k<=n+1;k+=1){
            ptf = 1.;
            for (i=1;i<=n;i+=1){
                ptf += dptf;
                sum = E_Ei(Aridge, ptf, etaf, phif)*dptf;
                
                // cout<<pti<<std::setw(7)<<yi<<std::setw(9)<<phii<<endl;       
                // cout<<mb*mb+pti*pti<<setw(10)<<sqrt(mb*mb+pti*pti)<<endl;
                
            }
            // cout<<etaf<<setw(15)<<phif<<setw(15)<<ptf<<setw(15)<<sum<<setw(15)<<lightcone(ptf,yi)<<endl;
            // fprintf(fpt, "%f, %f, %f\n", etaf, phif, sum);
            fout2<<etaf<<","<<phif<<","<<sum<<endl;
            phif += dphif;
            sum = 0.;
        }
        etaf += detaf;
    }

    fout.close();


    std::ofstream fout3("E_Ei/rapidity_eta=1.5.csv");
    fout3<<"ptf,yi,yf,cosh(yf)/cosh(yi)\n";

    etaf = 1.5;
    phif = 0.;
    ptf = 0.15;
    dptf = double ((4.-.15)/1000);

    for(i=1;i<=1000;i++){
        double pti = sqrt(ptf*ptf-2*ptf*q*cos(phif)/cosh(etajet)+q*q/(cosh(etajet)*cosh(etajet)));
        double yi = rapidityintit(pti,etaf);
        double yf = rapidityintit(ptf,etaf);
        fout3<<ptf<<','<<yi<<','<<yf<<','<<cosh(yf)/cosh(yi)<<endl;
        ptf += dptf;
    }

    fout3.close();
*/




/*

 //3차원 함수 연습  etaf, phif
    // cout<<'2'<<endl;

    double dptf, plots;
    fpt = fopen("3dplot_Ridge_Jet.csv", "w+");
    fprintf(fpt,"etaf, phif, dF\n");
    // double result[10][1000] = {0.};
    

    // double detaf, dphif, etaf, phif, ptf;

    n = 100;
    dptf = double ((4-0.15)/n);
    detaf = double ((5.+5.)/n);
    dphif = double ((1.+1.)/n);
    etaf = -5.;
    phif = -1.;
    double ptf0 = 0.15;
    sum = 0.;


    check2 = 1;
    double constant = Njet/(Tjet*(m+Tjet)*2*M_PI);

    for(j=1;j<=n+1;j+=1){
        phif = -1.;
        for(k=1;k<=n+1;k+=1){
            ptf = 0.15;
            for (i=1;i<=n;i+=1){
                ptf += dptf;
                double yi = rapidityintit(sqrt(ptf*ptf-2*ptf*q*cos(phif)+q*q),etaf);
                check = plots;
                plots = ptf*RidgeDis(Aridge, ptf, etaf, phif, check2)*fRNk/3+fj*integralNjet(ptf, etaf, phif, constant);
                // plots = ptf*RidgeDis(Aridge, ptf, etaf, phif, check2);
                // plots = RidgeDis(Aridge, ptf, etaf, phif, check2);

                sum += plots;
                

                if (i==n){
                    sum += ptf0*(RidgeDis(Aridge, ptf0, etaf, phif, check2)*fRNk/3+fj*integralNjet(ptf0, etaf, phif, constant)-plots)/2.;
                    // sum += (ptf0*RidgeDis(Aridge, ptf0, etaf, phif, check2)-plots)/2.;
                    // sum += (RidgeDis(Aridge, ptf0, etaf, phif, check2)-plots)/2.;
                    sum *= dptf;
                    break;
                }
                else if (i !=1 && lightcone(sqrt((ptf+dptf)*(ptf+dptf)-2*(ptf+dptf)*q*cos(phif)+q*q),yi)>=1.){
                    sum += (ptf0*RidgeDis(Aridge, ptf0, etaf, phif, check2)*fRNk/3+fj*integralNjet(ptf0, etaf, phif, constant)-plots)/2.;
                    // sum += (ptf0*RidgeDis(Aridge, ptf0, etaf, phif, check2)-plots)/2.;     //Ridge 계산할 경우
                    // sum += (RidgeDis(Aridge, ptf0, etaf, phif, check2)-plots)/2.;     //E/Ei 계산할 경우
                    sum *= dptf;
                    // cout<<"!!!!"<<endl;
                    break;
                }

                // cout<<etaf<<setw(15)<<phif<<setw(15)<<ptf<<setw(15)<<sum<<setw(15)<<lightcone(ptf,yi)<<endl;

                
                // cout<<pti<<std::setw(7)<<yi<<std::setw(9)<<phii<<endl;       
                // cout<<mb*mb+pti*pti<<setw(10)<<sqrt(mb*mb+pti*pti)<<endl;
                
            }
            // cout<<etaf<<setw(15)<<phif<<setw(15)<<ptf<<setw(15)<<sum<<setw(15)<<lightcone(ptf,yi)<<endl;
            fprintf(fpt, "%f, %f, %f\n", etaf, phif, sum);
            phif += dphif;
        }
        etaf += detaf;
    }


 //3차원 함수 연습 etaf, ptf




    // // double dptf, plots;
    // fpt = fopen("3d_etapt_E.csv", "w+");
    // fprintf(fpt,"etaf, ptf, dF\n");
    // // double result[10][1000] = {0.};
    // // int n;

    // // double detaf, dphif, etaf, phif, ptf, phif0;
    // // double phif0;

    // n = 100;
    // dptf = double ((4-0.15)/n);
    // detaf = double ((5.+5.)/n);
    // dphif = double ((1.+1.)/n);
    // etaf = -5.;
    // phif = -1.;
    // ptf0 = 0.15;
    // sum = 0.;
    // phif0 = -1.;




    // // double constant = Njet/(Tjet*(m+Tjet)*2*M_PI);

    // for(j=1;j<=n+1;j+=1){
    //     ptf = 0.15;
    //     for(k=1;k<=n+1;k+=1){
    //         phif = -1.;
    //         for (i=1;i<=n;i+=1){
    //             phif += dphif;
    //             double yi = rapidityintit(sqrt(ptf*ptf-2*ptf*q*cos(phif)+q*q),etaf);
    //             check = plots;
    //             // plots = ptf*RidgeDis(Aridge, ptf, etaf, phif, check2)*fRNk/3+fj*integralNjet(ptf, etaf, phif, constant);
    //             plots = ptf*RidgeDis(Aridge, ptf, etaf, phif, check2);
    //             // plots = RidgeDis(Aridge, ptf, etaf, phif, check2);

    //             sum += plots;
                

    //             if (i==n){
    //                 // sum += ptf*(RidgeDis(Aridge, ptf, etaf, phif, check2)*fRNk/3+fj*integralNjet(ptf, etaf, phif0, constant)-plots)/2.;
    //                 sum += (ptf*RidgeDis(Aridge, ptf, etaf, phif, check2)-plots)/2.;
    //                 // sum += (RidgeDis(Aridge, ptf, etaf, phif, check2)-plots)/2.;
    //                 sum *= dphif;
    //                 break;
    //             }
    //             else if (i !=1 && lightcone(ptf-1,yi)>=1.){
    //                 // sum += (ptf*R=np.loadtxt('Test.csv',delimiter=',',usecols=[0],skiprows=1)
    //                 sum *= dptf;
    //                 // cout<<"!!!!"<<endl;
    //                 break;
    //             }

    //             //lightcone은 phi에 관련이 없으므로 제외.

    //             // cout<<etaf<<setw(15)<<ptf<<setw(15)<<phif<<setw(30)<<sum<<setw(15)<<lightcone(ptf-1,yi)<<endl;

                
    //             // cout<<pti<<std::setw(7)<<yi<<std::setw(9)<<phii<<endl;       
    //             // cout<<mb*mb+pti*pti<<setw(10)<<sqrt(mb*mb+pti*pti)<<endl;
                
    //         }
    //         // cout<<etaf<<setw(15)<<phif<<setw(15)<<ptf<<setw(15)<<sum<<setw(15)<<lightcone(ptf,yi)<<endl;
    //         fprintf(fpt, "%f, %f, %f\n", etaf, ptf, sum);
    //         ptf += dptf;
    //     }
    //     etaf += detaf;
    // }

    // //eta에 대해서만 그려보자.

    // fpt = fopen("eta_E.csv", "w+");
    // totalsum = 0.;
    // n = 100;
    // dptf = double ((4-0.15)/n);
    // detaf = double ((5.+5.)/n);
    // dphif = double ((1.+1.)/n);
    // etaf = -5.;
    // phif = -1.;
    // ptf0 = 0.15;
    // sum = 0.;
    // phif0 = -1.;


    // for(j=1;j<=n+1;j+=1){
    //     phif = -1.;
    //     for(k=1;k<=n+1;k+=1){
    //         ptf = 0.15;
    //         for (i=1;i<=n;i+=1){
    //             ptf += dptf;
    //             double yi = rapidityintit(ptf-1,etaf);
    //             check = plots;
    //             // plots = ptf*RidgeDis(Aridge, ptf, etaf, phif, check2)*fRNk/3+fj*integralNjet(ptf, etaf, phif, constant);
    //             plots = ptf*RidgeDis(Aridge, ptf, etaf, phif, check2);
    //             // plots = RidgeDis(Aridge, ptf, etaf, phif, check2);

    //             sum += plots;
                

    //             if (i==n){
    //                 // sum += ptf0*(RidgeDis(Aridge, ptf0, etaf, phif, check2)*fRNk/3+fj*integralNjet(ptf0, etaf, phif, constant)-plots)/2.;
    //                 sum += (ptf0*RidgeDis(Aridge, ptf0, etaf, phif, check2)-plots)/2.;
    //                 // sum += (RidgeDis(Aridge, ptf0, etaf, phif, check2)-plots)/2.;
    //                 sum *= dptf;
    //                 break;
    //             }
    //             else if (i !=1 && lightcone(ptf+d    fpt = fopen("dNdy.csv", "w+");
    fprintf(fpt,"y, dFdy\n");

    n = 888;
    dpti = double ((8.8-0.0)/n);
    sum = 0.;


    pti = 0.;
    check2 = 0;

    for(yi=0.;yi<=5.4;yi+=0.1){
        pti = 0.;
        for(j=1;j<=n;j+=1){
            
            sum += Aridge*integralAridge(pti, yi, check2)*dpti;
            pti += dpti;
        }
        fprintf(fpt,"%f, %f\n", yi, sum);
        sum = 0.;
    }
                
    //             // cout<<pti<<std::setw(7)<<yi<<std::setw(9)<<phii<<endl;       
                
    //         }
    //         // cout<<etaf<<setw(15)<<phif<<setw(1    fpt = fopen("dNdy.csv", "w+");
    fprintf(fpt,"y, dFdy\n");

    n = 888;
    dpti = double ((8.8-0.0)/n);
    sum = 0.;


    pti = 0.;
    check2 = 0;

    for(yi=0.;yi<=5.4;yi+=0.1){
        pti = 0.;
        for(j=1;j<=n;j+=1){
            
            sum += Aridge*integralAridge(pti, yi, check2)*dpti;
            pti += dpti;
        }
        fprintf(fpt,"%f, %f\n", yi, sum);
        sum = 0.;
    }
    // cout<<ptf<<setw(15)<<sum<<setw(15)<<lightcone(ptf,yi)<<endl;
    //         totalsum += sum*dphif;
    //         phif += dphif;
    //     }
    //     fprintf(fpt, "%f, %f\n", etaf, totalsum);
    //     totalsum = 0.;
    //     etaf += detaf;
    // }




    //pt고정, 3dplot
    
    
    fpt = fopen("pt_3dplot_Ridge.csv", "w+");
    fprintf(fpt,"etaf, phif, pt=0.15, pt=1, pt=2, pt=3, pt=4\n");
    // double result[10][1000] = {0.};
    // int n;

    // double detaf, dphif, etaf, phif, ptf;

    n = 1000;
    // dptf = double ((4-0.15)/n);
    detaf = double ((1.4+1.4)/n);
    dphif = double ((2.5+2.5)/n);
    etaf = -1.4;
    phif = -2.5;
    // ptf = 4;
    
    double Ni = 1.;
    double sigmay = 5.5;
    double Ai = Ni*exp(m/T)/(pow(2*M_PI,3/2)*sigmay*T);

    check2 = 1;

    for(j=1;j<=n+1;j+=1){
        phif = -2.5;
        for(k=1;k<=n+1;k+=1){
            // cout<<etaf<<setw(15)<<phif<<setw(15)<<ptf<<setw(15)<<sum<<setw(15)<<lightcone(ptf,yi)<<endl;
            fprintf(fpt, "%f, %f, %f, %f, %f, %f, %f\n", etaf, phif, RidgeDis(Aridge, 0.15, etaf, phif, check2), RidgeDis(Aridge, 1, etaf, phif, check2), RidgeDis(Aridge, 2, etaf, phif, check2), RidgeDis(Aridge, 3, etaf, phif, check2), RidgeDis(Aridge, 4, etaf, phif, check2));
            phif += dphif;
        }
        etaf += detaf;
    }

*/


    fclose(fpt);



    return 0;
}
