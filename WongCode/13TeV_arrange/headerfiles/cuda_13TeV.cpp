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
#include "headerfiles/function.hpp"



void ptyield(double ptq){
    
    std::ofstream foutalicept("pTdis.csv");
    foutalicept<<"pt,ALICE_Ridge,CMS_Ridge,ALICE_Jet,CMS_Jet\n";

    double delta_Deltaeta = 0.4;
    double delta_Deltaetacms = 4.;

    std::cout<<"CZYAM Calculate start"<<std::endl;

    std::cout<<"CZYAM Calculate start"<<std::endl;


    // double CZYAM[2][1100]={0.};
    int n_ptyield = 400;
    double *CZYAM = (double*)malloc(2*n_ptyield*sizeof(double));
    for(int k=0;k<=n_ptyield;k++){
        for(int l=1;l<=n_1;l++){
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



}