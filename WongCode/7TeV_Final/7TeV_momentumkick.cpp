#define _USE_MATH_DEFINES
// #define _WIN32_WINNT 0x0A00

#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <fstream>
#include <thread>


//Ridge parameters
double a = .5;    //fall off parameter
double T = .63;    //Temperature, GeV
double q = .88;    //GeV
double m = 0.13957018;  //m == mpi
double mb = m; //mb==mpi, GeV
double md = 1.;   //Ge+V
double sqrSnn = 7000.;
double mp = 0.938272046; //Proton mass, GeV

double x = .6;
double y = .5;

// double A = 630.;
// double C = 6.084;
// double B = 0.9;

double fRNk01_1 = 2.0;

double fRNk1_2 = 2.;
double fRNk2_3 = 5.5;
double fRNk3_4 = 15.125;


double fRNk1_4 = .9;
double fRNk = 4.;// double x = 0.99

double frnkconst = 7.;

double etajet = 0.;



//Jet Parameters
double Njet=.75;
double fj=1.;
double Tjet=0.6;
double sigmaphizero = .5;
double ma = .5;
// double ma = 10000.;


double Aridge;
double constant = Njet/(Tjet*(m+Tjet)*2*M_PI);

double frnk(double pt){
    // return exp(-frnkconst*pt/(frnkconst*frnkconst+pt*pt));
    // return exp(-frnkconst*pt/sqrt(md*md+pt*pt));
    // return frnkconst*exp(-pt/sqrt(md*md+pt*pt));
    // return (md*exp(-pt/sqrt(md*md+pt*pt)))/frnkconst;
    // return exp(-md*pt/sqrt(md*md+pt*pt));
    // return  3.5*pt-3.25;
    return  x*exp(y*pt);
}

double rapidityintit(double pt, double eta){
    double root=sqrt((pt*pt*cosh(eta)*cosh(eta))+mb*mb);
    double a = root+pt*sinh(eta);
    double b = root-pt*sinh(eta);
    return (log(a/b))/2;
}

double lightcone(double pti, double yi){
    // double yi = rapidityinit(pti);
    double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mp
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

    // std::cout<<ptf<<std::setw(20)<<phif<<std::setw(20)<<etaf<<std::setw(20)<<integralAridge(pti, yi, check)<<std::endl;
    // std::cout<<pti<<std::setw(20)<<yi<<std::setw(20)<<x<<std::endl;

    if(check == 1){
        if(x>=1.){
            return 0.;
        }
        else{
            // std::cout<<pti<<std::setw(20)<<yi<<std::setw(20)<<pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti)<<std::endl;
            return pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
        }
        
    }

    //Aridge
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
    double ptisq = ptf*ptf-2*ptf*q*cos(phif)/cosh(etajet)+q*q/(cosh(etajet)*cosh(etajet));
    double pti = sqrt(ptisq);
    if(ptisq<0.0000000001){
        pti = 0.;
    }
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

    // std::cout<<ptf<<std::setw(20)<<phif<<std::setw(20)<<etaf<<std::setw(20)<<pti<<std::setw(20)<<Ei<<std::endl;
    if(ptf*ptf-2*ptf*q*cos(phif)+q*q<0){
        // std::cout<<std::endl<<std::endl;
    }

    if (x>=1.){
        return 0.;
    }
    
    else{              
        return (Aridge*integralAridge(pti, yi, check))*sqrt(1.-((mb*mb)/((mb*mb+ptf*ptf)*cosh(yf)*cosh(yf))))*(E/Ei);           // E/Ei 임을 명심하자.
    }//
    
}

//Jetdistribution
double integralNjet(double pt, double eta, double phi, double constant){
    double sigmaphi;
    sigmaphi = (sigmaphizero*sigmaphizero*ma*ma)/(ma*ma+pt*pt);
    return (constant/sigmaphi)*exp(((m-sqrt(m*m+pt*pt))/Tjet)-((phi*phi+eta*eta)/(2*sigmaphi)));
}
//여기까지 Jet


//phi correlation 적분 범위

int n_2 = 100;
int n_3 = 100;
int n_4 = 100;
int n_5 = 100;
int n_6 = 100;

double delta_Deltaeta = 0.4;
double delta_Deltaetacms = 4.;

double dptf_2 = double ((2.-1.)/n_2);   //func2 pt 적분 범위 : 1~2
double ptf0_2 = 1.;
double dptf_3 = double ((3.-2.)/n_3);   //func3 pt 적분 범위 : 2~3
double ptf0_3 = 2.;
double dptf_4 = double ((4.-3.)/n_4);   //func4 pt 적분 범위 : 3~4
double ptf0_4 = 3.;
double dptf_5 = double ((4.-1.)/n_5);   //func5 pt 적분 범위 : 1~4
double ptf0_5 = 1.;

double dptf_6 = double ((0.1-1.)/n_6);   //func6 pt 적분 범위 : 0.1~1
double ptf0_6 = .1;

//나중에 따로 하게될 수 있으니 귀찮더라도 하나씩 따로 분리하자.

double detaf_2 = double ((1.8-1.6)/n_2);  //func2 eta 적분 범위 : 1.6~1.8
double etaf0_2 = 1.6;
double detacms_2 = double ((4.-2.)/n_2);  //func2 eta 적분 범위 : 2~4(cms)
double etacms0_2 = 2.;

double detaf_3 = double ((1.8-1.6)/n_3);  //func3 eta 적분 범위 : 1.6~1.8
double etaf0_3 = 1.6;
double detacms_3 = double ((4.-2.)/n_2);  //func3 eta 적분 범위 : 2~4(cms)
double etacms0_3 = 2.;

double detaf_4 = double ((1.8-1.6)/n_4);  //func4 eta 적분 범위 : 1.6~1.8
double etaf0_4 = 1.6;
double detacms_4 = double ((4.-2.)/n_2);  //func4 eta 적분 범위 : 2~4(cms)
double etacms0_4 = 2.;

double detaf_5 = double ((1.8-1.6)/n_5);  //func5 eta 적분 범위 : 1.6~1.8
double etaf0_5 = 1.6;
double detacms_5 = double ((4.-2.)/n_2);  //func5 eta 적분 범위 : 2~4(cms)
double etacms0_5 = 2.;

double detaf_6 = double ((4.-2.)/n_6);  //func6 eta 적분 범위 : 2~4
double etaf0_6 = 2.;


double dphif_2 = double ((1.28+1.28)/n_2);  //func2 phi 출력 범위 : -1.28~1.28
double phif_2 = -1.28;
double dphif_3 = double ((1.28+1.28)/n_3);  //func3 phi 출력 범위 : -1.28~1.28
double phif_3 = -1.28;
double dphif_4 = double ((1.28+1.28)/n_4);  //func4 phi 출력 범위 : -1.28~1.28
double phif_4 = -1.28;
double dphif_5 = double ((1.28+1.28)/n_5);  //func5 phi 출력 범위 : -1.28~1.28
double phif_5 = -1.28;

double dphif_6 = double ((1.28+1.28)/n_6);  //func6 phi 출력 범위 : -1.28~1.28
double phif_6 = -1.28;

int check2_2 = 1;
int check2_3 = 1;
int check2_4 = 1;
int check2_5 = 1;

int check2_6 = 1;

//pt distribution 적분 범위

int n_1 = 50;
double detaf_1 = double((1.8-1.6)/n_1);     //func1 eta 적분범위(Ridge, alice) : 1.6~1.8 -> x2 해야 함.
double etaf0_1 = 1.6;
double detacms_1 = double((4.-2.)/n_1);     //func1 eta 적분범위(Ridge, cms) : 1.6~1.8 -> x2 해야 함.
double etacms0_1 = 2.;

double detajet_1 = double((1.6-0.)/n_1);     //func1 eta 적분범위(Jet) : 0~1.6 -> x2 해야 함.
double etajet0_1 = 0.;

double dphif_1 = double((1.28+1.28)/n_1);   //func1 phi 적분범위 : -1.28~1.28
double phif_1 = -1.28;
double phif0_1 = -1.28;
double delta_Deltaphi = 2.56;

int n_1_pt = 20;        //pt 포인트당 적분범위 이므로 주의하자. 현재 0.5간격으로 적분중
double dptf_1 = double ((10.-0.)/n_1_pt);  //func1 pt 출력 범위 : 0.15~11
double ddptf_1 = dptf_1/n_1;
double ptf_1 = 0.;
int check2_1 = 1;

void func1(){
    double etaf_1, dist_1, sum_1_alice, sum_1_cms, sum_1j, ptf_11, etajet_1, etacms_1, sum_1_alice_czyam, sum_1_cms_czyam;
    int i_1, j_1, k_1, l_1, h_1;
    
    std::ofstream foutalicept("pTdis.csv");
    foutalicept<<"pt,ALICE_Ridge,CMS_Ridge,ALICE_Jet,CMS_Jet\n";
    double CZYAM[30][2]={0.};

    // cout<<"1"<<endl;
    
    for(k_1=1;k_1<=n_1_pt+1;k_1++){
        for(l_1=1;l_1<=n_1;l_1++){
            etaf_1 = etaf0_1;
            etacms_1 = etacms0_1;
            for (i_1=1;i_1<=n_1+1;i_1++){
                sum_1_alice += ptf_1*frnk(ptf_1)*RidgeDis(Aridge, ptf_1, etaf_1, 1.28, check2_1)*detaf_1*ddptf_1/delta_Deltaeta;
                sum_1_cms += ptf_1*frnk(ptf_1)*RidgeDis(Aridge, ptf_1, etacms_1, 1.28, check2_1)*detacms_1*ddptf_1/delta_Deltaetacms;
                etaf_1 += detaf_1;
                etacms_1 += detacms_1;
            }
            ptf_1 += ddptf_1;
        }
        sum_1_alice *= 2*2/3;
        sum_1_cms *= 2*2/3;
        CZYAM[k_1][0]=sum_1_alice;
        CZYAM[k_1][1]=sum_1_cms;
        // std::cout<<ptf_1<<std::setw(20)<<sum_1_alice<<std::endl;
        sum_1_alice = sum_1_cms = 0.;
    }


    sum_1_alice = sum_1_cms = sum_1_alice_czyam = sum_1_alice_czyam = sum_1j = 0.;



    for(k_1=1;k_1<=n_1_pt;k_1++){
        phif_1 = phif0_1;
        for(j_1=1;j_1<=n_1+1;j_1++){
            sum_1_alice_czyam += CZYAM[k_1][0]*dphif_1;
            sum_1_cms_czyam += CZYAM[k_1][1]*dphif_1;
            phif_1 += dphif_1;
        }
        CZYAM[k_1][0] = sum_1_alice_czyam;
        CZYAM[k_1][1] = sum_1_cms_czyam;
        sum_1_alice_czyam = sum_1_cms_czyam = 0.;
    }


    double ptf_1 = 0.;
    double phif_1 = phif0_1;

    for(k_1=1;k_1<=n_1_pt+1;k_1++){
        for(l_1=1;l_1<=n_1;l_1++){
            // std::cout<<ptf_1<<std::endl;
            phif_1 = phif0_1;
            // std::cout<<ptf-q<<std::endl;
            for(j_1=1;j_1<=n_1+1;j_1++){
                etaf_1 = etaf0_1;
                etacms_1 = etacms0_1;
                etajet_1 = etajet0_1;
                // for (i_1=1;i_1<=n_1+1;i_1++){
                for (i_1=1;i_1<=n_1+1;i_1++){
                    
                    // etaf += detaf;
                    // double yi = rapidityintit(sqrt(ptf*ptf-2*ptf*q*cos(phif)+q*q),etaf+detaf);
                    // check = dist;
                    // dist = RidgeDis(Aridge, ptf, etaf, phif, check2);

                    // // cout<<etaf<<endl;
                    
                    // sum += dist;
                    // dist_1 = RidgeDis(Aridge, ptf_1, etaf_1, phif_1, check2_1)*detaf_1*dphif_1/delta_Deltaeta;   //임시저장 용도
                    // sum_1 += RidgeDis(Aridge, ptf_1, etaf_1, phif_1, check2_1)*detaf_1*dphif_1/delta_Deltaeta;
                    double dist = ptf_1*frnk(ptf_1)*RidgeDis(Aridge, ptf_1, etaf_1, phif_1, check2_1)*detaf_1*dphif_1*ddptf_1/delta_Deltaeta;
                    sum_1_alice += dist;
                    // sum_1_alice += 
                    sum_1_cms += ptf_1*frnk(ptf_1)*RidgeDis(Aridge, ptf_1, etacms_1, phif_1, check2_1)*detacms_1*dphif_1*ddptf_1/delta_Deltaetacms;
                    sum_1j += ptf_1*integralNjet(ptf_1,etajet_1,phif_1,constant)*detajet_1*dphif_1*ddptf_1/delta_Deltaphi;
                    // std::cout<<ptf_1<<std::setw(20)<<phif_1<<std::setw(20)<<etaf_1<<std::setw(20)<<RidgeDis(Aridge, ptf_1, etaf_1, phif_1, check2_1)<<std::endl;
                    // std::cout<<ptf_1<<std::setw(20)<<phif_1<<std::setw(20)<<etaf_1<<std::setw(20)<<sum_1_alice<<std::endl;
                    // std::cout<<ptf_1<<std::setw(5);
                    etaf_1 += detaf_1;
                    etacms_1 += detacms_1;
                    etajet_1 += detajet_1;
                }

                // totalsum += sum*dphif;
                phif_1 += dphif_1;
                // std::cout<<ptf_1<<std::setw(20)<<phif_1<<std::setw(20)<<sum_1<<std::endl;
                // sum = 0.;
            }
            // std::cout<<ptf_1<<std::setw(20)<<sum_1<<std::endl;
            // sum_1 = (sum_1*2*fRNk)/(3*2*M_PI);
            // sum_1 = 2*(sum_1*2*fRNk)/3;

            ptf_1 += ddptf_1;
        }
        // ptf_1 -= dptf_1;
        // std::cout<<std::endl<<ptf_1<<std::endl;

        // std::cout<<ptf_1<<std::endl; sum_1_alice_czyam = sum_1_cms_czyam =
        sum_1_alice *= 2*2/3;
        sum_1_cms *= 2*2/3;
        sum_1j *= 2;
        
        // std::cout<<sum_1_alice_czyam<<std::endl;
        // double ptfout = ptf_1-(dptf_1/2);
        foutalicept<<ptf_1-(dptf_1/2)<<","<<sum_1_alice-CZYAM[k_1][0]<<','<<sum_1_cms<<","<<sum_1j<<std::endl;
        // std::cout<<ptf_1<<std::endl;

        //-CZYAM[k_1][1]

        // sum_1_alice = sum_1_cms = sum_1j = 0.;
        sum_1_alice = sum_1_cms = sum_1_alice_czyam = sum_1_cms_czyam = sum_1j = 0.;
        // ptf_1 += dptf_1;
        // std::cout<<ptf_1<<std::endl;
        // resultsum = resultsum + totalsum;
        // cout<<pt<<std::setw(20)<<totalsum<<endl;

    }
    foutalicept.close();
}


//pt1~2, 2~3, 3~4 부터 시작해보자.


void func2(){
    
    double ptf_2, etaf_2, sum_2_alice, sum_2_cms, sum_2j, etacms_2;
    int i_2, j_2, k_2;
    std::ofstream fout12("phiCorrelation_pt1-2.csv");
    fout12<<"phi,Alice,CMS\n";
    
    sum_2_alice = sum_2_cms = sum_2j = 0.;
    
    for(k_2=1;k_2<=n_2+1;k_2++){
        etacms_2 = etacms0_2;
        etaf_2 = etaf0_2;
        // cout<<phif<<endl;
        for(i_2=1;i_2<=n_2+1;i_2++){
            
            ptf_2 = ptf0_2;
            // cout<<etaf<<endl;
            for(j_2=1;j_2<=n_2+1;j_2+=1){
                
                sum_2_alice += frnk(ptf_2)*ptf_2*RidgeDis(Aridge, ptf_2, etaf_2, phif_2, check2_2)*dptf_2*detaf_2/delta_Deltaeta;
                sum_2_cms += frnk(ptf_2)*ptf_2*RidgeDis(Aridge, ptf_2, etacms_2, phif_2, check2_2)*dptf_2*detacms_2/delta_Deltaetacms;
                
                // sum_2 += ptf_2*RidgeDis(Aridge, ptf_2, etaf_2, phif_2, check2_2)*dptf_2*detaf_2/delta_Deltaeta;
                // sum_2j += ptf_2*integralNjet(ptf_2,etaf_2,phif_2,constant)*dptf_2*detaf_2;
                // sum_2 *= frnk(ptf_2);
                // std::cout<<phif_2<<std::setw(15)<<etaf_2<<std::setw(15)<<ptf_2<<std::setw(15)<<sum_2<<std::endl;
                ptf_2 += dptf_2;
            }
            
            etaf_2 += detaf_2;
            etacms_2 += detacms_2;
            // std::cout<<i_2<<std::setw(10)<<etaf_2<<std::endl;
        }
        // sum_2 *= 2*2*fRNk1_2/3;
        sum_2_alice *= 2*2/3;
        sum_2_cms *= 2*2/3;
        // sum_2j *= fj;
        // fprintf(fpt,"%f, %f\n", phif, sum);
        

        fout12<<phif_2<<","<<sum_2_alice<<","<<sum_2_cms<<std::endl;
        phif_2 += dphif_2;
        sum_2_alice = sum_2_cms = sum_2j = 0.;
    }
    fout12.close();
}

void func3(){
    double ptf_3, etaf_3, sum_3_alice, sum_3_cms, sum_3j, etacms_3;
    int i_3, j_3, k_3;
    std::ofstream fout23 ("phiCorrelation_pt2-3.csv");
    fout23<<"phi,ALICE,CMS\n";

    sum_3_alice = sum_3_cms = sum_3j = 0.;

    for(k_3=1;k_3<=n_3+1;k_3++){
        etacms_3 = etacms0_3;
        etaf_3 = etaf0_3;
        // cout<<phif<<endl;
        for(i_3=1;i_3<=n_3+1;i_3++){
            ptf_3 = ptf0_3;
            // cout<<etaf<<endl;
            for(j_3=1;j_3<=n_3+1;j_3+=1){
                sum_3_alice += frnk(ptf_3)*ptf_3*RidgeDis(Aridge, ptf_3, etaf_3, phif_3, check2_3)*dptf_3*detaf_3/delta_Deltaeta;
                sum_3_cms += frnk(ptf_3)*ptf_3*RidgeDis(Aridge, ptf_3, etacms_3, phif_3, check2_3)*dptf_3*detacms_3/delta_Deltaetacms;
                // sum_3 += ptf_3*RidgeDis(Aridge, ptf_3, etaf_3, phif_3, check2_3)*dptf_3*detaf_3/delta_Deltaeta;
                sum_3j += ptf_3*integralNjet(ptf_3,etaf_3,phif_3,constant)*dptf_3*detaf_3;
                
                // cout<<phif<<setw(15)<<etaf<<setw(15)<<ptf<<setw(15)<<sum<<endl;
                ptf_3 += dptf_3;
            }
            etaf_3 += detaf_3;
            etacms_3 += detacms_3;
        }
        // sum_3 *= 2*2*fRNk2_3/3;
        sum_3_alice *= 2*2/3;
        sum_3_cms *= 2*2/3;
        sum_3j *= fj;
        // fprintf(fpt,"%f, %f\n", phif, sum);

        fout23<<phif_3<<","<<sum_3_alice<<","<<sum_3_cms<<std::endl;
        sum_3_alice = sum_3_cms = sum_3j = 0.;
        phif_3 += dphif_3;
    }
    fout23.close();
}

void func4(){
    double ptf_4, etaf_4, sum_4_alice, sum_4_cms, sum_4j, etacms_4;
    int i_4, j_4, k_4;
    std::ofstream fout34("phiCorrelation_pt3-4.csv");
    fout34<<"phi,ALICE,CMS\n";

    sum_4_alice = sum_4_cms = sum_4j = 0.;

    for(k_4=1;k_4<=n_4+1;k_4++){
        etacms_4 = etacms0_4;
        etaf_4 = etaf0_4;
        // cout<<phif<<endl;
        for(i_4=1;i_4<=n_4+1;i_4++){
            ptf_4 = ptf0_4;
            // cout<<etaf<<endl;
            for(j_4=1;j_4<=n_4+1;j_4+=1){
                sum_4_alice += frnk(ptf_4)*ptf_4*RidgeDis(Aridge, ptf_4, etaf_4, phif_4, check2_4)*dptf_4*detaf_4/delta_Deltaeta;
                sum_4_cms += frnk(ptf_4)*ptf_4*RidgeDis(Aridge, ptf_4, etacms_4, phif_4, check2_4)*dptf_4*detacms_4/delta_Deltaetacms;
                // sum_4 += ptf_4*RidgeDis(Aridge, ptf_4, etaf_4, phif_4, check2_4)*dptf_4*detaf_4/delta_Deltaeta;
                // sum_4j += ptf_4*integralNjet(ptf_4,etaf_4,phif_4,constant)*dptf_4*detaf_4;
                // sum_4 *= frnk(ptf_4);
                // cout<<phif<<setw(15)<<etaf<<setw(15)<<ptf<<setw(15)<<sum<<endl;
                ptf_4 += dptf_4;
            }
            etaf_4 += detaf_4;
            etacms_4 += detacms_4;
        }
        // sum_4 *= 2*2*fRNk3_4/3;
        sum_4_alice *= 2*2/3;
        sum_4_cms *= 2*2/3;
        sum_4j *= fj;
        // fprintf(fpt,"%f, %f\n", phif, sum);

        fout34<<phif_4<<","<<sum_4_alice<<","<<sum_4_cms<<std::endl;
        sum_4_alice = sum_4_cms = 0.;
        phif_4 += dphif_4;
    }
    fout34.close();
}


void func5(){
    double ptf_5, etaf_5, sum_5_alice, sum_5_cms, sum_5j, etacms_5;
    int i_5, j_5, k_5;
    std::ofstream fout14("phiCorrelation_pt1-4.csv");
    fout14<<"phi,ALICE,CMS\n";

    sum_5_alice = sum_5_cms = sum_5j = 0.;

    for(k_5=1;k_5<=n_5+1;k_5++){
        etaf_5 = etaf0_5;
        etacms_5 = etacms0_5;
        // cout<<phif<<endl;
        for(i_5=1;i_5<=n_5+1;i_5++){
            ptf_5 = ptf0_5;
            // cout<<etaf<<endl;
            for(j_5=1;j_5<=n_5+1;j_5+=1){
                // sum_5 += frnk(ptf_5)*ptf_5*RidgeDis(Aridge, ptf_5, etaf_5, phif_5, check2_5)*dptf_5*detaf_5;
                sum_5_alice += frnk(ptf_5)*ptf_5*RidgeDis(Aridge, ptf_5, etaf_5, phif_5, check2_5)*dptf_5*detaf_5/delta_Deltaeta;
                sum_5_cms += frnk(ptf_5)*ptf_5*RidgeDis(Aridge, ptf_5, etacms_5, phif_5, check2_5)*dptf_5*detacms_5/delta_Deltaetacms;
                // sum_5j += ptf_5*integralNjet(ptf_5,etaf_5,phif_5,constant)*dptf_5*detaf_5;
                // sum_5 *= frnk(ptf_5);
                // cout<<phif<<setw(15)<<etaf<<setw(15)<<ptf<<setw(15)<<sum<<endl;
                ptf_5 += dptf_5;
            }
            etaf_5 += detaf_5;
            etacms_5 += detacms_5;
        }
        // sum_5 *= 2*2*fRNk1_4/3;
        sum_5_alice *= 2*2/3;
        sum_5_cms *= 2*2/3;
        sum_5j *= fj;
        // fprintf(fpt,"%f, %f\n", phif, sum);

        fout14<<phif_5<<","<<sum_5_alice<<","<<sum_5_cms<<std::endl;
        sum_5_alice = sum_5_cms = sum_5j = 0.;
        phif_5 += dphif_5;
    }
    fout14.close();        
}

void func6(){
    double ptf_6, etaf_6, sum_6, sum_6j;
    int i_6, j_6, k_6;
    std::ofstream fout011("phiCorrelation_pt01-1.csv");
    fout011<<"phi,Ridge,Jet\n";

    sum_6 = sum_6j = 0.;


    for(k_6=1;k_6<=n_6;k_6++){
        etaf_6 = etaf0_6;
        // cout<<phif<<endl;
        for(i_6=1;i_6<=n_6;i_6++){
            ptf_6 = ptf0_6;
            // cout<<etaf<<endl;
            for(j_6=1;j_6<=n_6;j_6+=1){
                // sum_6 += frnk(ptf_6)*ptf_6*RidgeDis(Aridge, ptf_6, etaf_6, phif_6, check2_6)*dptf_6*detaf_6;
                sum_6 += ptf_6*RidgeDis(Aridge, ptf_6, etaf_6, phif_6, check2_6)*dptf_6*detaf_6;
                sum_6j += ptf_6*integralNjet(ptf_6,etaf_6,phif_6,constant)*dptf_6*detaf_6;
                // cout<<phif<<setw(15)<<etaf<<setw(15)<<ptf<<setw(15)<<sum<<endl;
                ptf_6 += dptf_6;
            }
            etaf_6 += detaf_6;
        }
        sum_6 *= 2*2*fRNk01_1/3;
        // sum_6 *= 2*2/3;
        sum_6j *= fj;
        // fprintf(fpt,"%f, %f\n", phif, sum);

        fout011<<phif_6<<","<<sum_6<<","<<sum_6j<<std::endl;
        sum_6 = sum_6j = 0.;
        phif_6 += dphif_6;
    }
    fout011.close();
}

int n_7 = 100;

double dptf_7 = double ((40.-10.)/n_7);
double ptf0_7 = 10.;
double detaf_7 = double ((1.8-1.6)/n_7);
double etaf0_7 = 1.6;
double dphif_7 = double ((1.28+1.28)/n_7);  //func2 phi 출력 범위 : -1.28~1.28
double phif_7 = -1.28;
int check2_7 = 1;






void func7(){
    double ptf_7, etaf_7, phif0_7, sum_7j, sum_7;
    int i_7, j_7, k_7;
    std::ofstream fout0002("Jet_PhiCorrelation_pT_over10.csv");
    fout0002<<"phi,d^2N/dphi\n";
    // double dptf;


    for(k_7=1;k_7<=n_7+1;k_7++){
        etaf_7 = etaf0_7;
        // cout<<phif<<endl;
        for(i_7=1;i_7<=n_7+1;i_7++){
            ptf_7 = ptf0_7;
            // cout<<etaf<<endl;
            for(j_7=1;j_7<=n_7+1;j_7++){
                // sum_7 = ptf_7*RidgeDis(Aridge, ptf_7, etaf_7, phif_7, check2_7)*1/3+fj*integralNjet(ptf_7, etaf_7, phif_7, constant);
                // sum_7 = frnk(ptf_7)*ptf_7*fj*integralNjet(ptf_7, etaf_7, phif_7, constant)*dphif_7*detaf_7;
                // sum_7 *= frnk(ptf_7);
                // sum_7 += ptf_7*RidgeDis(Aridge, ptf_7, etaf_7, phif_7, check2_7)*dphif_7*detaf_7;
                // cout<<ptf_7<<setw(15)<<etaf_7<<setw(15)<<phif_7<<setw(15)<<sum_7<<endl;
                sum_7j += ptf_7*integralNjet(ptf_7,etaf_7,phif_7,constant)*dptf_7*detaf_7/delta_Deltaeta;
                ptf_7 += dptf_7;
                
            }
            
            etaf_7 += detaf_7;
        }
        sum_7 *= 2*2/3;
        sum_7j *= 2;
        // fprintf(fpt,"%f, %f\n", phif, sum);

        fout0002<<phif_7<<","<<sum_7j<<std::endl;
        sum_7 = sum_7j = 0.;
        phif_7 += dphif_7;
    }
    fout0002.close();
}


void func8(){
    double ptf_8, etaf_8, phif_8, detaf_8, dphif_8, etaf0_8, phif0_8, dist_8, dptf_8, sum_8, sum_8j;
    int netaf_8, nphif_8, n_8, check2_8, i_8, j_8;
    std::ofstream fout1012("Y^near.csv");
    fout1012<<"pt^Jet,Y^near\n";
    // double dptf;

    n_8 = 100;
    // dptf_7 = double ((1.0-0.1)/n_7);
    detaf_8 = double ((1.2-1.0)/n_8);
    dphif_8 = double ((M_PI+M_PI)/n_8);
    sum_8 = 0.;
    etaf_8 = 1.;
    phif_8 = -M_PI;

    check2_8 = 1;

    for(ptf_8=0.1;ptf_8<=2.;ptf_8+=0.01){
        etaf_8 = 1.;
        // cout<<phif<<endl;
        for(i_8=1;i_8<=n_8+1;i_8++){
            phif_8 = -M_PI;
            // cout<<etaf<<endl;
            for(j_8=1;j_8<=n_8+1;j_8++){
                // sum_8 = ptf_8*RidgeDis(Aridge, ptf_8, etaf_8, phif_8, check2_8)*1/3+fj*integralNjet(ptf_8, etaf_8, phif_8, constant);
                // sum_8 = frnk(ptf_8)*ptf_8*fj*integralNjet(ptf_8, etaf_8, phif_8, constant)*dphif_8*detaf_8;
                sum_8 += ptf_8*RidgeDis(Aridge, ptf_8, etaf_8, phif_8, check2_8)*dphif_8*detaf_8;
                sum_8j += ptf_8*integralNjet(ptf_8,etaf_8,phif_8,constant)*dphif_8*dptf_8;
                // cout<<ptf_8<<setw(15)<<etaf_8<<setw(15)<<phif_8<<setw(15)<<sum_8<<endl;
                phif_8 += dphif_8;
            }
            etaf_8 += detaf_8;
        }
        sum_8 *= 2*2/3;
        // fprintf(fpt,"%f, %f\n", phif, sum);

        fout1012<<ptf_8<<","<<sum_8<<","<<sum_8j<<std::endl;
        sum_8 = sum_8j = 0.;
    }
    fout1012.close();
}

double n_9 = 200;
double detaf_9 = double((2.4+2.4)/n_9);     //func1 eta 적분범위(Ridge, alice) : 1.6~1.8 -> x2 해야 함.
double etaf0_9 = -2.4;

double dphif_9 = double((M_PI+M_PI)/n_9);   //func1 phi 적분범위 : -1.28~1.28
double phif_9 = -M_PI;
double phif0_9 = -M_PI;

double dptf_9 = double ((6.-0.1)/n_9);
double ptf_9 = 0.1;
int check2_9 = 1;

void func9(){
    double etaf_9, sum_9, sum_9j;
    int i_9, j_9, k_9;

    std::ofstream foutptdis("dNdetaptdpt.csv");
    foutptdis<<"pt,dN/detadpt_Ridge,dN/detadpt_Jet\n";

    sum_9 = sum_9j = 0.;

    for(k_9=1;k_9<=n_9+1;k_9++){
        phif_9 = phif0_9;
        for(j_9=1;j_9<=n_9+1;j_9++){
            etaf_9 = etaf0_9;
            for (i_9=1;i_9<=n_9+1;i_9++){
                
                // etaf += detaf;
                // double yi = rapidityintit(sqrt(ptf*ptf-2*ptf*q*cos(phif)+q*q),etaf+detaf);
                // check = dist;
                // dist = RidgeDis(Aridge, ptf, etaf, phif, check2);

                // // cout<<etaf<<endl;
                
                // sum += dist;
                sum_9 += frnk(ptf_9)*RidgeDis(Aridge, ptf_9, etaf_9, phif_9, check2_9)*detaf_9*dphif_9;
                sum_9j += integralNjet(ptf_9,etaf_9,phif_9,constant)*detaf_9*dphif_9;
                etaf_9 += detaf_9;
            }
            
            // totalsum += sum*dphif;
            phif_9 += dphif_9;

            
            // sum = 0.;

        }
        // cout<<ptf<<setw(20)<<sum<<endl;

        // sum_1 = (sum_1*2*fRNk)/(3*2*M_PI);
        sum_9 *= 2/(3*2*M_PI);
        sum_9j *= 1/(2*M_PI);

        foutptdis<<ptf_9<<","<<sum_9<<","<<sum_9j<<std::endl;

        ptf_9 += dptf_9;
        sum_9 = sum_9j = 0.;
        // std::cout<<ptf_1<<std::endl;
        // resultsum = resultsum + totalsum;
        // cout<<pt<<std::setw(20)<<totalsum<<endl;

    }
    foutptdis.close();
}


int main()
{
    //multithread 컴파일 할 경우 다음의 명령어를 입력하자.
    // g++ -o Ridge Ridge.cpp -lpthread
    using std::cout;
    using std::endl;
    using std::setw;
    using std::thread;
    
    // std::string buffer;
    double dyi, dphii, sum, totalsum, phii, yi, dpti, pti, sum2, resultsum;
    int i, j, k, nyi, npti, nphii, check2;

    // FILE *fpt;
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
    Aridge = 1/totalsum;
    // cout<<totalsum<<setw(15)<<Aridge<<endl;
    
    // double Aridge = 1/resultsum;
    cout<<totalsum<<setw(20)<<Aridge<<endl;
    
    //Aridge=0.1420386444451514

    // Aridge=0.1420386444451514;
    //여기까지 Aridge를 구하기 위한 적분

    

    //여기부터 Ridge항 적분 pt : 0.15~4, etaf : -1.4~+1.4, phif : -1~1 ->func1

    std::ofstream fout("FittingParameters_Ridge.csv");

    fout<<"a,"<<a<<endl;
    fout<<"T,"<<T<<endl;
    fout<<"q,"<<q<<endl;
    fout<<"pion mass,"<<m<<endl;
    fout<<"beam mass,"<<mb<<endl;
    fout<<"m_d,"<<md<<endl;
    fout<<"sqrtSnn,"<<sqrSnn<<endl;
    fout<<"A_Ridge,"<<Aridge<<endl;
    // fout<<"k,"<<frnkconst<<",f_R<N_k>=k/sqrt(k^2+pt)"<<endl;
    // fout<<"k,"<<frnkconst<<",f_R<N_k>=k*exp(-pt/sqrt(md*md+pt*pt))"<<endl;
    // fout<<"f_R<N_k>(pt=0.1~1),"<<fRNk01_1<<endl;
    // fout<<"f_R<N_k>(pt=1~2),"<<fRNk1_2<<endl;
    // fout<<"f_R<N_k>(pt=2~3),"<<fRNk2_3<<endl;
    // fout<<"f_R<N_k>(pt=3~4),"<<fRNk3_4<<endl;
    // fout<<"f_R<N_k>(pt=1~4),"<<fRNk1_4<<endl;
    // fout<<"f_R<N_k>=3.5pt-3.25,"<<endl;
    fout<<"f_R<N_k>=,"<<x<<"exp("<<y<<"*pT)"<<endl;


    // std::ofstream foutjet("FittingParameters_Jet.csv");

    // foutjet<<"N_jet,"<<Njet<<endl;
    // foutjet<<"f_j,"<<0.632<<endl;
    // foutjet<<"T_jet,"<<Tjet<<endl;
    // foutjet<<"sigma_phizero,"<<sigmaphizero<<endl;
    // foutjet<<"m_a,"<<ma<<endl;


    fout.close();
    
    // foutjet.close();
    
    // func1();
    // func2();
    // func3();
    // func4();
    // func5();
    
    // func6();

    thread t1(func1);    //pt distribution - alice
    // thread t9(func9);    //pt distribution - cms

    thread t2(func2);   //1D eta correlation, pt = 1~2
    thread t3(func3);   //1D eta correlation, pt = 2~3
    thread t4(func4);   //1D eta correlation, pt = 3~4
    thread t5(func5);   //1D eta correlation, pt = 1~4
    thread t9(func9);   //dN/pTdpT   -> 포기
    // thread t7(func7);
    // thread t6(func6);
    
    // thread t8(func8);
    
    t1.join();
    t9.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    // t9.join();
    // t7.join();

    // t6.join();
    
    // t8.join();

    char ch[100];
    std::string command = "python3 Graphs.py";
    strcpy(ch,command.c_str());
    std::system(ch);


    return 0;
}
