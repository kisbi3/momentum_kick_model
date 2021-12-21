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
#include "headerfiles/cuda_function.cuh"
#include "headerfiles/cuda_integral.cuh"
#include "headerfiles/integral.hpp"
// #include "headerfiles/function.hpp"

int n_1 = 200;
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

int n_1_pt = 400;
double dptf_1 = double ((10.-0.)/n_1_pt);  //func1 pt 출력 범위 : 0.15~11
double ddptf_1 = dptf_1/n_1;
double ptf_1 = 0.;
int check2_1 = 1;

void ptyield(double ptq){
    
    std::ofstream foutalicept("pTdis.csv");
    foutalicept<<"pt,ALICE_Ridge,CMS_Ridge,ALICE_Jet,CMS_Jet\n";

    double delta_Deltaeta = 0.4;
    double delta_Deltaetacms = 4.;

    // 'sum_1_alice', 'sum_1_cms' in origin code
    double ptczyam_alice, ptczyam_cms;
    ptczyam_alice = ptczyam_cms = 0.;

    std::cout<<"CZYAM Calculate start"<<std::endl;

    // double CZYAM[2][1100]={0.};

    //pt_yield number of points('n_1_pt' in origin code)
    int n_ptyield = 400;
    double ptf_end = 10.;
    double ptf_st = 0.;
    double dptf = (ptf_end-ptf_st)/double(n_ptyield);
    double ptf_distend = ptf_end/double(n_ptyield);

    //pt_yield number of integral bin ('n_1' in origin code)
    int n_ptyield_int = 200;

    double eta_alice_st = 1.6;
    double eta_alice_end = 1.8;
    double eta_cms_st = 2.;
    double eta_cms_end = 4.;

    double *CZYAM = (double*)malloc(2*n_ptyield*sizeof(double));
    for(int k=0;k<=n_ptyield;k++){
        CZYAM[0][k] = 2*(2/3)*(cuda_second(ptf_st+k*dptf, ptf_distend+k*dptf, eta_alice_st, eta_elice_end, 2)/delta_Deltaeta);
        CZYAM[1][k] = 2*(2/3)*(cuda_second(ptf_st+k*dptf, ptf_distend+k*dptf, eta_cms_st, eta_cms_end, 2)/delta_Deltaetacms);
        // ptczyam_alice = cuda_second(ptf_st+k*dptf, ptf_distend+k*dptf, eta_alice_st, eta_elice_end, 2)/delta_Deltaeta;
        // ptczyam_cms = cuda_second(ptf_st+k*dptf, ptf_distend+k*dptf, eta_cms_st, eta_cms_end, 2)/delta_Deltaetacms;

        // ptczyam_alice *= 2*2/3;
        // ptczyam_cms *= 2*2/3;
        // CZYAM[0][k]=ptczyam_alice;
        // CZYAM[1][k]=ptczyam_cms;
        ptczyam_alice = ptczyam_cms = 0.;
    }

    std::cout<<"CZYAM Calculate 1/2..."<<std::endl;

    double phif = -1.28;
    double dist1, dist2;
    dist1 = dist2 = 0.;
    for(int k=0;k<=n_ptyield;k++){
        phif = -1.28;
        for(int j=0;j<=n_ptyield_int;j++){
            dis1 += CZYAM[0][k]*dphif_1;
            dist2 += CZYAM[1][k]*dphif_1;
            phif += double((1.28+1.28)/n_ptyield_int);
        }
        CZYAM[0][k] = dis1;
        CZYAM[1][k] = dist2;
        dist1 = distt2 = 0.;
    }


    std::cout<<"CZYAM Calculate end"<<std::endl;


    for(int k=1;k<=(n_1_pt+2)/2;k_1++){
        std::thread t11(func11, ptf_1, phif_1, k_1, ptq,'r', ddptf_1);
        std::thread t12(func11, ptf_2, phif_1, k_2, ptq,'r', ddptf_1);
        t11.join();
        t12.join();
        ptf_1 += dptf_1;
        ptf_2 += dptf_1;
        // std::cout<<Yridge[0][k_1-1]<<std::setw(20)<<Yridge[1][k_1-1]<<std::endl;
        k++;
    }

    norm[0] = 1./norm[0];
    norm[1] = 1./norm[1];
    // std::cout<<norm[0]<<std::setw(20)<<norm[1]<<std::endl;
    ptf_1=0.;
    for(k_1=1;k_1<=n_1_pt+1;k_1++){
        // std::cout<<k_1<<std::endl;
        // std::cout<<k_1<<std::setw(20)<<Yridge[0][k_1-1]<<std::setw(20)<<Yridge[1][k_1-1]<<std::endl;
        foutalicept<<ptf_1+(dptf_1/2)<<","<<Yridge[0][k_1-1]*norm[0]<<","<<Yridge[1][k_1-1]*norm[1]<<","<<sum_1j<<std::endl;
        ptf_1 += dptf_1;
    }
    foutalicept.close();

}