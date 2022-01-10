#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <fstream>
#include <thread>
#include <time.h>
#include <mutex>


//Ridge parameters
double a = .5;    //fall off parameter
double T = .65;    //Temperature, GeV
double q_1 = .9;    //GeV
double m = 0.13957018;  //m == mpi
double mb = m; //mb==mpi, GeV
double md = 1.;   //GeV
double sqrSnn = 13000.;
double mp = 0.938272046; //Proton mass, GeV


double xx = .67;
double yy = .71;


double etajet = 0.;


//Jet Parameters
double Njet=10.;
double fj=1.;
double Tjet=0.7;
double sigmaphizero = .05;
double ma = 100.;
// double ma = 10000.;


double Aridge;
double constant = Njet/(Tjet*(m+Tjet)*2*M_PI);

double frnk(double pt){

    return  xx*exp(yy*pt);

    // return exp(-frnkconst*pt/(frnkconst*frnkconst+pt*pt));
    // return exp(-frnkconst*pt/sqrt(md*md+pt*pt));
    // return frnkconst*exp(-pt/sqrt(md*md+pt*pt));
    // return (md*exp(-pt/sqrt(md*md+pt*pt)))/frnkconst;
    // return exp(-md*pt/sqrt(md*md+pt*pt));
    // return  3.5*pt-3.25;

    // if(0.<=pt && pt<1.){
    //     return xx*exp(yy*0.5);
    // }
    // else if(1.<=pt && pt<2.){
    //     return xx*exp(yy*1.5);
    // }
    // else if(2.<=pt && pt<3.){
    //     return xx*exp(yy*2.5);
    // }
    // else if(3.<=pt && pt<4.){
    //     return xx*exp(yy*3.5);
    // }
    // else if(4.<=pt && pt<5.){
    //     return xx*exp(yy*4.5);
    // }
    // else if(5.<=pt && pt<6.){
    //     return xx*exp(yy*5.5);
    // }
    // else if(6.<=pt && pt<7.){
    //     return xx*exp(yy*6.5);
    // }
    // else if(7.<=pt && pt<8.){
    //     return xx*exp(yy*7.5);
    // }
    // else if(8.<=pt && pt<9.){
    //     return xx*exp(yy*8.5);
    // }
    // else if(9.<=pt && pt<=11.){
    //     return xx*exp(yy*10.);
    // }
    // else{
    //     printf("\nSomething error, ptf : %f\n", pt);
    //     exit(1);
    //     return 0;
    // }


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
double RidgeDis(double Aridge, double ptf, double etaf, double phif, int check, double q1){
    double ptisq = ptf*ptf-2*ptf*q1*cos(phif)/cosh(etajet)+q1*q1/(cosh(etajet)*cosh(etajet));
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
    if(ptf*ptf-2*ptf*q1*cos(phif)+q1*q1<0){
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



//pt distribution 적분 범위

int n_1 = 100;
double detaf_1 = double((1.8-1.6)/n_1);     //func1 eta 적분범위(Ridge, alice) : 1.6~1.8 -> x2 해야 함.
double etaf0_1 = 1.6;
double lasteta = 1.8;
double detacms_1 = double((4.-2.)/n_1);     //func1 eta 적분범위(Ridge, cms) : 1.6~1.8 -> x2 해야 함.
double etacms0_1 = 2.;
double lastetacms = 4.;

double detajet_1 = double((1.6-0.)/n_1);     //func1 eta 적분범위(Jet) : 0~1.6 -> x2 해야 함.
double etajet0_1 = 0.;
double lastetajet = 1.6;

double dphif_1 = double((1.28+1.28)/n_1);   //func1 phi 적분범위 : -1.28~1.28
double phif_1 = -1.28;
double phif0_1 = -1.28;
double lastphi = 1.28; //적분 끝 범위
double delta_Deltaphi = 2.56;

int n_1_pt = 200;
double dptf_1 = double ((15.-0.15)/n_1_pt);  //func1 pt 출력 범위 : 0.15~11
double ddptf_1 = dptf_1/n_1;
double ptf_1 = 0.15;
int check2_1 = 1;

// double arr[3]= {0.};

//multithread 검색했을 때에 void가 아닌 경우가 없었음 -> return 안되는듯?
//이번 경우에는 계속 더하기만 하면 되기 때문에 mutex로 더할 필요 없음. 단, 다른 메모리에서 접근하여 가져갈 위험이 있는 경우, mutex의 lock/unlock 필수!

void func10(double arr[], double ptf, double phif, double etaf, double etacms, double etajet, double qq, double dptff){
    int i;
    double delta_Deltaeta = 0.4;
    double delta_Deltaetacms = 4.;
    for(i=1;i<=n_1/2;i++){
        arr[0] += ptf*frnk(ptf)*RidgeDis(Aridge, ptf, etaf, phif, check2_1, qq)*detaf_1*dphif_1*dptff/delta_Deltaeta;
        arr[1] += ptf*frnk(ptf)*RidgeDis(Aridge, ptf, etacms, phif, check2_1, qq)*detacms_1*dphif_1*dptff/delta_Deltaetacms;
        arr[2] += ptf*integralNjet(ptf,etajet,phif,constant)*detajet_1*dphif_1*dptff/delta_Deltaphi;
        etaf += detaf_1;
        etacms += detacms_1;
        etajet += detajet_1;
    }
}

void func9(double arr[], double ptf_1, double phif_1, double q, double dptfff){
    int i_9, j_9;
    double etaf_1, etacms_1, etajet_1;
    for(j_9=1;j_9<=n_1/2;j_9++){
        etaf_1 = etaf0_1;
        etacms_1 = etacms0_1;
        etajet_1 = etajet0_1;
        std::thread t11(func10,arr, ptf_1, phif_1, etaf_1, etacms_1, etajet, q, dptfff);
        std::thread t12(func10,arr, ptf_1, phif_1, (etaf_1+lasteta)/2., (etacms_1+lastetacms)/2., (etajet_1+lastetajet)/2., q, dptfff);
        // std::thread t13(func10,arr, ptf_1, phif_1, (etaf_1+lasteta)/2., (etacms_1+lastetacms)/2., (etajet_1+lastetajet)/2., q, dptfff);
        t11.join();
        t12.join();
        // t13.join();
        phif_1 += dphif_1;
    }
    // std::cout<<arr[0]<<std::setw(15)<<arr[1]<<std::endl;
}

double CZYAM[2][1100]={0.};
double Yridge[2][1100]={0.};
double CZYAMj[2][1100]={0.};
double Yridgej[2][1100]={0.};
double norm[2] = {0.};
double normj[2] = {0.};

void func11(double ptf2, double phif, int k, double q, char dist, double dptf2){
    // std::cout<<dist<<std::endl;
    int l;
    double arr1[3] = {0.};
    double arr2[3] = {0.};
    // for(l=1;l<=n_1;l++){
    //     double phif_1 = phif0_1;
    //     double phif_2 = (phif0_1+lastphi)/2.;
    //     // phif_1 = phif0_1;
    //     std::thread t9(func9,arr, ptf2, phif_1, q);
    //     std::thread t10(func9,arr, ptf2, phif_2, q);
    //     t9.join();
    //     t10.join();
    //     ptf2 += ddptf_1;
    // }
    // std::cout<<ptf2<<std::setw(20)<<arr[0]<<std::setw(20)<<arr[1]<<std::endl;
    // // std::cout<<ptf_1<<std::endl;
    // arr[0] *= 2.*2./3.;
    // arr[1] *= 2.*2./3.;
    // arr[2] *= 2.*2./3.;

    if(dist == 'r'){
        for(l=1;l<=n_1;l++){
            double phif_r1 = phif0_1;
            double phif_r2 = (phif0_1+lastphi)/2.;
            // phif_1 = phif0_1;
            std::thread t9(func9,arr1, ptf2, phif_r1, q, dptf2);
            std::thread t10(func9,arr1, ptf2, phif_r2, q, dptf2);
            //ddptf_1
            t9.join();
            t10.join();
            ptf2 += dptf2;
        }


        std::cout<<"Ridge"<<std::setw(20)<<ptf2<<std::setw(20)<<arr1[0]<<std::setw(20)<<arr1[1]<<std::endl;


        // std::cout<<ptf_1<<std::endl;
        arr1[0] *= 2.*2./3.;
        arr1[1] *= 2.*2./3.;
        arr1[2] *= 2.*2./3.;
        Yridge[0][k-1] = arr1[0]-CZYAM[0][k];
        Yridge[1][k-1] = arr1[1];


        // std::cout<<k-1<<std::setw(20)<<Yridge[0][k-1]<<std::setw(20)<<Yridge[1][k-1]<<std::endl;

        //Alice의 경우 1<pt<4의 경우에서 그래프가 그려진다.
        if(1.<=ptf2 && ptf2<=4.){
            norm[0] += (arr1[0]-CZYAM[0][k])*dptf_1;   
            // norm[1] += arr1[1]*dptf_1;                     
        }
        // norm[0] += (arr1[0]-CZYAM[0][k])*dptf_1;   
        // norm[1] += (arr1[1]-CZYAM[1][k])*dptf_1;   //cms도 alice와 동일하게 해야되는거 아닌가? 왜이러지?

        //CMS의 경우에는 0.154~14.37임.
        if(0.154<=ptf2 && ptf2<=14.37){
            norm[1] += arr1[1]*dptf_1;    
        }
        
    }

    else if(dist == 'j'){
        for(l=1;l<=n_1;l++){
            double phif_j1 = phif0_1;
            double phif_j2 = (phif0_1+lastphi)/2.;
            // phif_1 = phif0_1;
            std::thread t14(func9,arr2, ptf2, phif_j1, q, dptf2);
            std::thread t15(func9,arr2, ptf2, phif_j2, q, dptf2);
            t14.join();
            t15.join();
            // std::cout<<ptf2<<std::endl;
            ptf2 += dptf2;
        }


        
        std::cout<<"Jet"<<std::setw(20)<<q<<std::setw(20)<<arr2[0]<<std::setw(20)<<arr2[1]<<std::endl;
        
        
        // std::cout<<ptf_1<<std::endl;
        arr2[0] *= 2.*2./3.;
        arr2[1] *= 2.*2./3.;
        arr2[2] *= 2.*2./3.;
        Yridgej[0][k-1] = arr2[0]-CZYAMj[0][k];
        // Yridgej[0][k-1] = arr2[0];
        Yridgej[1][k-1] = arr2[1];
        //std::cout<<k-1<<std::setw(20)<<Yridgej[0][k-1]<<std::setw(20)<<Yridgej[1][k-1]<<std::endl;

        // if(1.<=ptf2 && ptf2<=4.){
        //     normj[0] += (arr2[0]-CZYAMj[0][k])*dptf_1;   
        //     normj[1] += arr2[1]*dptf_1;                     
        // }        

        normj[0] += (arr2[0]-CZYAMj[0][k])*dptf_1;   
        normj[1] += arr2[1]*dptf_1;                     
    }

}



void func1(double ptq){
    double etaf_1, dist_1, sum_1_alice, sum_1_cms, sum_1j, ptf_11, etajet_1, etacms_1, sum_1_alice_czyam, sum_1_cms_czyam;
    int i_1, j_1, k_1, l_1, h_1;
    // double** Yridge = new double*[n_1_pt];    //normalization 하자. 동적배열 delete 필수!
    // double** CZYAM = new double*[n_1_pt];
    // for(int i=0;i<2;i++){
    //     Yridge[i] = new double[n_1_pt];
    //     CZYAM[i] = new double[n_1_pt];
    // }

    // double norm[2] = {0.};    //normalization

    
    
    std::ofstream foutalicept("pTdis.csv");
    foutalicept<<"pt,ALICE_Ridge,CMS_Ridge,ALICE_Jet,CMS_Jet\n";
    // double CZYAM[2][400]={0.};
    // double Yridge[2][400]={0.};

    double delta_Deltaeta = 0.4;
    double delta_Deltaetacms = 4.;
    
    std::cout<<"CZYAM Calculate start"<<std::endl;

    for(k_1=1;k_1<=n_1_pt+1;k_1++){
        for(l_1=1;l_1<=n_1;l_1++){
            etaf_1 = etaf0_1;
            etacms_1 = etacms0_1;
            for (i_1=1;i_1<=n_1+1;i_1++){
                sum_1_alice += ptf_1*frnk(ptf_1)*RidgeDis(Aridge, ptf_1, etaf_1, 1.28, check2_1, ptq)*detaf_1*ddptf_1/delta_Deltaeta;
                sum_1_cms += ptf_1*frnk(ptf_1)*RidgeDis(Aridge, ptf_1, etacms_1, 1.28, check2_1, ptq)*detacms_1*ddptf_1/delta_Deltaetacms;
                etaf_1 += detaf_1;
                etacms_1 += detacms_1;
            }
            ptf_1 += ddptf_1;
        }
        sum_1_alice *= 2*2/3;
        sum_1_cms *= 2*2/3;
        CZYAM[0][k_1]=sum_1_alice;
        CZYAM[1][k_1]=sum_1_cms;
        sum_1_alice = sum_1_cms = 0.;
    }
    //phi Correlation 처럼 일단 pt, eta적분을 통해(phi=1.28) 가장 작은 값 계산

    std::cout<<"CZYAM Calculate 1/2..."<<std::endl;

    sum_1_alice = sum_1_cms = sum_1_alice_czyam = sum_1_alice_czyam = sum_1j = 0.;


    for(k_1=1;k_1<=n_1_pt;k_1++){
        phif_1 = phif0_1;
        for(j_1=1;j_1<=n_1+1;j_1++){
            sum_1_alice_czyam += CZYAM[0][k_1]*dphif_1;
            sum_1_cms_czyam += CZYAM[1][k_1]*dphif_1;
            phif_1 += dphif_1;
        }
        CZYAM[0][k_1] = sum_1_alice_czyam;
        CZYAM[1][k_1] = sum_1_cms_czyam;
        sum_1_alice_czyam = sum_1_cms_czyam = 0.;
    }
    //phi에 대해서 적분
    std::cout<<"CZYAM Calculate end"<<std::endl;

    double ptf_1 = 0.15;
    // double ptf_2 = dptf_1*n_1_pt/2;
    double ptf_2 = (15.-0.15)/2.;
    double phif_1 = phif0_1;
    // double arr[3]= {0.};
    int k_2 = (n_1_pt+2)/2;

    for(k_1=1;k_1<=(n_1_pt+2)/2;k_1++){

        std::thread t11(func11, ptf_1, phif_1, k_1, ptq,'r', ddptf_1);
        std::thread t12(func11, ptf_2, phif_1, k_2, ptq,'r', ddptf_1);
        t11.join();
        t12.join();
        ptf_1 += dptf_1;
        ptf_2 += dptf_1;
        // std::cout<<Yridge[0][k_1-1]<<std::setw(20)<<Yridge[1][k_1-1]<<std::endl;
        k_2++;
    }

    // norm[0] = 1./norm[0];
    // norm[1] = 1./norm[1];


    // std::cout<<norm[0]<<std::setw(20)<<norm[1]<<std::endl;
    ptf_1=0.15;
    for(k_1=1;k_1<=n_1_pt+1;k_1++){
        // std::cout<<k_1<<std::endl;
        // std::cout<<k_1<<std::setw(20)<<Yridge[0][k_1-1]<<std::setw(20)<<Yridge[1][k_1-1]<<std::endl;

        // foutalicept<<ptf_1+(dptf_1/2)<<","<<Yridge[0][k_1-1]*norm[0]<<","<<Yridge[1][k_1-1]*norm[1]<<","<<sum_1j<<std::endl;

        //normalize 하면 안된다. 이때 계산한 값들의 적분값과 데이터의 적분값이 동일하도록 하자. 이는 graph 그리는 곳에서 진행한다.
        foutalicept<<ptf_1+(dptf_1/2)<<","<<Yridge[0][k_1-1]<<","<<Yridge[1][k_1-1]<<","<<sum_1j<<std::endl;
        ptf_1 += dptf_1;
    }
    foutalicept.close();

    std::ofstream fout("pTdis_integral.csv");
    fout<<"Alice,"<<norm[0]<<","<<"CMS,"<<norm[1]<<std::endl;
    fout.close();

    // delete [] Yridge;
}


//pt1~2, 2~3, 3~4 부터 시작해보자.

void func2(double ptf_st, double ptf_end, double etaf_st, double etaf_end, double etacms_st, double etacms_end, double etaatlas_st, double etaatlas_end, double phif_st, double phif_end, int n, int check2, double ptjetcut, double q){
    
    double ptf, etaf, phif, sum_alice, sum_cms, sum_j, etacms, etaatlas, sum_atlas, norm;
    double dptf, dphif, detaf, detacms, detaatlas, delta_Deltaeta, delta_Deltaetacms, delta_Deltaetaatlas, ptf0, etaf0, phif0, etacms0, etaatlas0;
    int i, j, k;
    norm = 0.;
    if(ptjetcut == 0.){
        std::string filename, pt_st, pt_end, filename2;
        filename = "phiCorrelation_pt";
        pt_st = std::to_string(int(ptf_st));
        pt_end = std::to_string(int(ptf_end));
        filename.append(pt_st);
        filename += "-";
        filename.append(pt_end);
        filename2 = filename;
        std::cout<<filename2<<std::setw(7)<<"start"<<std::endl;
        filename.append(".csv");

        std::ofstream fout(filename);
        fout<<"phi,Alice,CMS,ATLAS\n";

        dptf = (ptf_end-ptf_st)/n;
        detaf = (etaf_end-etaf_st)/n;
        detacms = (etacms_end-etacms_st)/n;
        detaatlas = (etaatlas_end-etaatlas_st)/n;
        dphif = (phif_end-phif_st)/n;

        ptf0 = ptf_st;
        etaf0 = etaf_st;
        etacms0 = etacms_st;
        etaatlas0 = etaatlas_st;
        phif0 = phif_st;
        
        sum_alice = sum_cms = sum_j = 0.;
        phif = phif0;
        
        delta_Deltaeta = 2*(etaf_end-etaf_st);
        delta_Deltaetacms = 2*(etacms_end-etacms_st);
        delta_Deltaetaatlas = 2*(etaatlas_end-etaatlas_st);


        for(k=1;k<=n+1;k++){
            etacms = etacms0;
            etaf = etaf0;
            etaatlas = etaatlas0;
            for(i=1;i<=n+1;i++){
                ptf = ptf0;
                for(j=1;j<=n+1;j+=1){
                    sum_alice += frnk(ptf)*ptf*RidgeDis(Aridge, ptf, etaf, phif, check2, q)*dptf*detaf/delta_Deltaeta;
                    sum_cms += frnk(ptf)*ptf*RidgeDis(Aridge, ptf, etacms, phif, check2, q)*dptf*detacms/delta_Deltaetacms;
                    sum_atlas += frnk(ptf)*ptf*RidgeDis(Aridge, ptf, etaatlas, phif, check2, q)*dptf*detaatlas/delta_Deltaetaatlas;
                    // std::cout<<ptf<<std::setw(20)<<FRNK<<std::endl;
                    ptf += dptf;
                }
                // std::cout<<etaatlas<<std::endl;
                
                etaf += detaf;
                etacms += detacms;
                etaatlas += detaatlas;
            }
            sum_alice *= 2*2/3;
            sum_cms *= 2*2/3;
            sum_atlas *= 2*2/3;

            norm += sum_atlas*dphif;
            
            fout<<phif<<","<<sum_alice<<","<<sum_cms<<","<<sum_atlas<<std::endl;
            phif += dphif;
            sum_alice = sum_cms = sum_atlas = sum_j = 0.;
        }
        if(phif <= 1.1){
            std::cout<<"phif : "<<phif<<"\tIntegrate result\t"<<norm<<std::endl;
        }
        fout.close();

        std::cout<<filename2<<std::setw(7)<<"end"<<std::endl;
    }
    else{
        //
        // q = q*ptjetcut;
        // q=1.;
        //
        q = 0.05*ptjetcut+0.5;
        std::string filename, ptcut, filename2;
        filename = "phiCorrelation_pt_jetcut_";
        ptcut = std::to_string(int(ptjetcut));
        filename.append(ptcut);
        filename2 = filename;
        std::cout<<filename2<<std::setw(7)<<"start"<<std::endl;
        filename.append(".csv");

        std::ofstream fout(filename);
        fout<<"phi,Alice,CMS\n";

        dptf = (ptf_end-ptf_st)/n;
        detaf = (etaf_end-etaf_st)/n;
        detacms = (etacms_end-etacms_st)/n;
        dphif = (phif_end-phif_st)/n;

        ptf0 = ptf_st;
        etaf0 = etaf_st;
        etacms0 = etacms_st;
        phif0 = phif_st;
        
        sum_alice = sum_cms = sum_j = 0.;
        phif = phif0;
        
        delta_Deltaeta = 2*(etaf_end-etaf_st);
        delta_Deltaetacms = 2*(etacms_end-etacms_st);


        for(k=1;k<=n+1;k++){
            etacms = etacms0;
            etaf = etaf0;
            for(i=1;i<=n+1;i++){
                ptf = ptf0;
                for(j=1;j<=n+1;j+=1){
                    
                    sum_alice += frnk(ptf)*ptf*RidgeDis(Aridge, ptf, etaf, phif, check2, q)*dptf*detaf/delta_Deltaeta;
                    sum_cms += frnk(ptf)*ptf*RidgeDis(Aridge, ptf, etacms, phif, check2, q)*dptf*detacms/delta_Deltaetacms;
                    ptf += dptf;
                }
                
                etaf += detaf;
                etacms += detacms;
            }
            sum_alice *= 2*2/3;
            sum_cms *= 2*2/3;
            
            fout<<phif<<","<<sum_alice<<","<<sum_cms<<std::endl;
            phif += dphif;
            sum_alice = sum_cms = sum_j = 0.;
        }
        fout.close();

        std::cout<<filename2<<std::setw(7)<<"end"<<std::endl;
    }

}

void func3(double ptq, double ptj_st, double ptj_end){
    double etaf_1, dist_1, sum_1_alice, sum_1_cms, sum_1j, ptf_11, etajet_1, etacms_1, sum_1_alice_czyam, sum_1_cms_czyam, etaf, phif;
    int i_1, j_1, k_1, l_1, h_1;
    double ptj = ptj_st;
    double dptj = (ptj_end-ptj_st)/n_1_pt;
    // double** Yridge = new double*[n_1_pt];    //normalization 하자. 동적배열 delete 필수!
    // double** CZYAM = new double*[n_1_pt];
    // for(int i=0;i<2;i++){
    //     Yridge[i] = new double[n_1_pt];
    //     CZYAM[i] = new double[n_1_pt];
    // }

    // double norm[2] = {0.};    //normalization

    //jet에 대해서 event cut하기. 그래도 그냥 Y_Ridge보단 쉬울듯? 그냥 ptjet에 대해서 적분하면 될듯. (1<p_T(assoc)<2)
    
    
    std::ofstream foutalicept("pTdis_jetcut.csv");
    foutalicept<<"pt,ALICE_Ridge,CMS_Ridge,ALICE_Jet,CMS_Jet\n";
    // double CZYAM[2][400]={0.};
    // double Yridge[2][400]={0.};

    double delta_Deltaeta = 0.4;
    // double delta_Deltaetacms = 4.;
    
    std::cout<<"CZYAM_Jet Calculate start"<<std::endl;

    double ptf, ptf0;
    ptf = ptf0 = 1.;
    double dptf = (2.-1.)/n_1;
    double detaf = (1.8-1.6)/n_1;
    double etaf0 = 1.6;
    int check2 = 1;
    double dphif = (1.28+1.28)/n_1;
    double phif0 = -1.28;


    for(k_1=1;k_1<=n_1_pt+1;k_1++){
        ptq = 0.05*ptj+0.5;
        ptf = ptf0;
        for(l_1=1;l_1<=n_1;l_1++){
            etaf = etaf0;
            // etacms_1 = etacms0_1;
            for (i_1=1;i_1<=n_1+1;i_1++){
                sum_1_alice += ptf*frnk(ptf)*RidgeDis(Aridge, ptf, etaf, 1.28, check2, ptq)*detaf*dptf/delta_Deltaeta;
                // sum_1_cms += ptf*frnk(ptf)*RidgeDis(Aridge, ptf, etacms_1, 1.28, check2_1, ptq)*detacms_1*ddptf_1/delta_Deltaetacms;
                etaf += detaf;
                // etacms_1 += detacms_1;
            }
            ptf += dptf;
        }
        sum_1_alice *= 2*2/3;
        // sum_1_cms *= 2*2/3;

        // std::cout<<ptj<<std::setw(20)<<sum_1_alice<<std::endl;

        CZYAMj[0][k_1]=sum_1_alice;
        // CZYAM[1][k_1]=sum_1_cms;
        sum_1_alice = sum_1_cms = 0.;
        ptj += dptj;
    }

    std::cout<<"CZYAM_Jet Calculate 1/2..."<<std::endl;

    sum_1_alice = sum_1_cms = sum_1_alice_czyam = sum_1_alice_czyam = sum_1j = 0.;


    for(k_1=1;k_1<=n_1_pt;k_1++){
        phif = phif0;
        for(j_1=1;j_1<=n_1+1;j_1++){
            sum_1_alice_czyam += CZYAMj[0][k_1]*dphif;
            // sum_1_cms_czyam += CZYAM[1][k_1]*dphif;
            phif += dphif;
        }
        CZYAMj[0][k_1] = sum_1_alice_czyam;
        // CZYAM[1][k_1] = sum_1_cms_czyam;
        sum_1_alice_czyam = sum_1_cms_czyam = 0.;
    }
    std::cout<<"CZYAM_Jet Calculate end"<<std::endl;

    // double ptf_1 = 0.;
    // double ptf_2 = dptf_1*n_1_pt/2;
    phif = phif0;
    
    ptj = ptj_st;
    double ptj2 = (ptj_end+ptj_st)/2.;
    // double arr[3]= {0.};
    int k_2 = (n_1_pt+2)/2;

    for(k_1=1;k_1<=(n_1_pt+2)/2;k_1++){
    // for(k_1=1;k_1<=n_1_pt+1;k_1++){
        ptq = 0.05*ptj+0.5;
        double ptq2 = 0.05*ptj2+0.5;
        ptf = ptf0;
        phif = phif0;
        std::thread t13(func11, ptf, phif, k_1, ptq, 'j', dptf);
        // std::cout<<ptj<<std::endl;
        std::thread t14(func11, ptf, phif, k_2, ptq2, 'j', dptf);
        t13.join();
        t14.join();
        ptj += dptj;
        ptj2 += dptj;
        // ptf_2 += dptf_1;
        // std::cout<<Yridge[0][k_1-1]<<std::setw(20)<<Yridge[1][k_1-1]<<std::endl;
        k_2++;
    }

    // norm[0] = 1./norm[0];
    // norm[1] = 1./norm[1];
    // std::cout<<norm[0]<<std::setw(20)<<norm[1]<<std::endl;
    ptj=ptj_st;
    for(k_1=1;k_1<=n_1_pt;k_1++){
        // std::cout<<k_1<<std::endl;
        // std::cout<<k_1<<std::setw(20)<<Yridge[0][k_1-1]<<std::setw(20)<<Yridge[1][k_1-1]<<std::endl;
        // foutalicept<<ptj<<","<<Yridge[0][k_1-1]*norm[0]<<","<<Yridge[1][k_1-1]*norm[1]<<","<<sum_1j<<std::endl;
        foutalicept<<ptj<<","<<Yridgej[0][k_1-1]<<","<<Yridgej[1][k_1-1]<<","<<sum_1j<<std::endl;
        ptj += dptj;
    }
    foutalicept.close();
    // delete [] Yridge;
}

int main()
{
    //multithread 컴파일 할 경우 다음의 명령어를 입력하자.
    // g++ -o Ridge Ridge.cpp -lpthread
    using std::cout;
    using std::endl;
    using std::setw;
    using std::thread;
    time_t start, end;
    // clock_t start, end;

    start = time(NULL);
    // start = clock();



    // std::string buffer;
    double dyi, dphii, sum, totalsum, phii, yi, dpti, pti, sum2, resultsum;
    int i, j, k, nyi, npti, nphii, check2;

    
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
    fout<<"q,"<<q_1<<endl;
    fout<<"pion mass,"<<m<<endl;
    fout<<"beam mass,"<<mb<<endl;
    fout<<"m_d,"<<md<<endl;
    fout<<"sqrtSnn,"<<sqrSnn<<endl;
    fout<<"A_Ridge,"<<Aridge<<endl;
    fout<<"f_R<N_k>=,"<<xx<<"exp("<<yy<<"*pT)"<<endl;


    // std::ofstream foutjet("FittingParameters_Jet.csv");

    // foutjet<<"N_jet,"<<Njet<<endl;
    // foutjet<<"f_j,"<<0.632<<endl;
    // foutjet<<"T_jet,"<<Tjet<<endl;
    // foutjet<<"sigma_phizero,"<<sigmaphizero<<endl;
    // foutjet<<"m_a,"<<ma<<endl;


    fout.close();

    // foutoh.close();
    
    // foutjet.close();
    

    thread t1(func1, q_1);   //pt distribution - alice
    thread t9(func3, q_1, 0., 50.);



    thread t2(func2, 1., 2., 1.6, 1.8, 2., 4., 2., 5., -1.28, 1.28, 300, 1, 0., q_1);   //1D phi correlation, pt = 1~2
    thread t3(func2, 2., 3., 1.6, 1.8, 2., 4., 2., 5., -1.28, 1.28, 300, 1, 0., q_1);   //1D phi correlation, pt = 2~3
    thread t4(func2, 3., 4., 1.6, 1.8, 2., 4., 2., 5., -1.28, 1.28, 300, 1, 0., q_1);   //1D phi correlation, pt = 3~4
    thread t5(func2, 1., 4., 1.6, 1.8, 2., 4., 2., 5., -1.28, 1.28, 300, 1, 0., q_1);   //1D phi correlation, pt = 1~4
    // q에 multiplicity dependence가 없는 경우
    thread t8(func2, .5, 5., 1.6, 1.8, 2., 4., 2., 5., -1., 1., 300, 1, 0., q_1);   //1D phi correlation, pt = 0.5~5
    
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t8.join();

    thread t6(func2, 1., 2., 1.6, 1.8, 2., 4., 2., 5., -1.28, 1.28, 300, 1, 10., q_1);   //1D phi correlation, pt = 1~2 (jet event cut)
    thread t7(func2, 1., 2., 1.6, 1.8, 2., 4., 2., 5., -1.28, 1.28, 300, 1, 20., q_1);   //1D phi correlation, pt = 1~2 (jet event cut)

    t7.join();
    t6.join();
    
    t1.join();
    t9.join();


    char ch[100];
    std::string command = "python3 Graphs.py";
    strcpy(ch,command.c_str());
    std::system(ch);

    // end = clock();
    end = time(NULL);
    cout<<double(end-start)<<endl;

    return 0;
}
