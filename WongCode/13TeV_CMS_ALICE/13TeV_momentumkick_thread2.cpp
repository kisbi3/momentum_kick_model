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


//Ridge parameters
double a = .5;    //fall off parameter
double T = .65;    //Temperature, GeV
double q = .9;    //GeV
double m = 0.13957018;  //m == mpi
double mb = m; //mb==mpi, GeV
double md = 1.;   //GeV
double sqrSnn = 13000.;
double mp = 0.938272046; //Proton mass, GeV


double x = .66;
double y = .71;


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


//pt distribution 적분 범위

int n_1 = 50;
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

int n_1_pt = 200;        //pt 포인트당 적분범위 이므로 주의하자. 현재 0.5간격으로 적분중
double dptf_1 = double ((10.-0.)/n_1_pt);  //func1 pt 출력 범위 : 0.15~11
double ddptf_1 = dptf_1/n_1;
double ptf_1 = 0.;
int check2_1 = 1;

// double arr[3]= {0.};

//multithread 검색했을 때에 void가 아닌 경우가 없었음 -> return 안되는듯?
//이번 경우에는 계속 더하기만 하면 되기 때문에 mutex로 더할 필요 없음. 단, 다른 메모리에서 접근하여 가져갈 위험이 있는 경우, mutex의 lock/unlock 필수!

void func10(double arr[], double ptf, double phif, double etaf, double etacms, double etajet){
    int i;
    double delta_Deltaeta = 0.4;
    double delta_Deltaetacms = 4.;
    for(i=1;i<=n_1/2;i++){
        arr[0] += ptf*frnk(ptf)*RidgeDis(Aridge, ptf, etaf, phif, check2_1)*detaf_1*dphif_1*ddptf_1/delta_Deltaeta;
        arr[1] += ptf*frnk(ptf)*RidgeDis(Aridge, ptf, etacms, phif, check2_1)*detacms_1*dphif_1*ddptf_1/delta_Deltaetacms;
        arr[2] += ptf*integralNjet(ptf,etajet,phif,constant)*detajet_1*dphif_1*ddptf_1/delta_Deltaphi;
        etaf += detaf_1;
        etacms += detacms_1;
        etajet += detajet_1;
    }
}

void func9(double arr[], double ptf_1, double phif_1){
    int i_9, j_9;
    double etaf_1, etacms_1, etajet_1;
    for(j_9=1;j_9<=n_1/2;j_9++){
        etaf_1 = etaf0_1;
        etacms_1 = etacms0_1;
        etajet_1 = etajet0_1;
        std::thread t11(func10,arr, ptf_1, phif_1, etaf_1, etacms_1, etajet);
        std::thread t12(func10,arr, ptf_1, phif_1, (etaf_1+lasteta)/2., (etacms_1+lastetacms)/2., (etajet_1+lastetajet)/2.);
        t11.join();
        t12.join();
        phif_1 += dphif_1;
    }
    // std::cout<<arr[0]<<std::setw(15)<<arr[1]<<std::endl;
}



void func1(){
    double etaf_1, dist_1, sum_1_alice, sum_1_cms, sum_1j, ptf_11, etajet_1, etacms_1, sum_1_alice_czyam, sum_1_cms_czyam;
    int i_1, j_1, k_1, l_1, h_1;
    double** Yridge = new double*[n_1_pt];    //normalization 하자. 동적배열 delete 필수!
    for(int i=0;i<2;i++){
        Yridge[i] = new double[n_1_pt];
    }
    double norm1 = 0., norm2 = 0.;    //normalization
    
    
    std::ofstream foutalicept("pTdis.csv");
    foutalicept<<"pt,ALICE_Ridge,CMS_Ridge,ALICE_Jet,CMS_Jet\n";
    double CZYAM[2][400]={0.};

    double delta_Deltaeta = 0.4;
    double delta_Deltaetacms = 4.;
    
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
        CZYAM[0][k_1]=sum_1_alice;
        CZYAM[1][k_1]=sum_1_cms;
        sum_1_alice = sum_1_cms = 0.;
    }

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


    double ptf_1 = 0.;
    double phif_1 = phif0_1;
    double arr[3]= {0.};

    for(k_1=1;k_1<=n_1_pt;k_1++){
        for(l_1=1;l_1<=n_1;l_1++){
            double phif_1_1 = phif0_1;
            double phif_1_2 = (phif0_1+lastphi)/2.;
            phif_1 = phif0_1;
            std::thread t9(func9,arr, ptf_1, phif_1_1);
            std::thread t10(func9,arr, ptf_1, phif_1_2);
            t9.join();
            t10.join();
            ptf_1 += ddptf_1;
        }
        std::cout<<ptf_1<<std::setw(20)<<arr[0]<<std::setw(20)<<arr[1]<<std::endl;
        // std::cout<<ptf_1<<std::endl;
        arr[0] *= 2.*2./3.;
        arr[1] *= 2.*2./3.;
        arr[2] *= 2.*2./3.;

        Yridge[0][k_1-1] = arr[0]-CZYAM[0][k_1];
        Yridge[1][k_1-1] = arr[1];

        if(1.<=ptf_1 && ptf_1<=4.){
            norm1 += (arr[0]-CZYAM[0][k_1])*dptf_1;   
            norm2 += arr[1]*dptf_1;                     
        }
        arr[0] = 0.;
        arr[1] = 0.;
        arr[2] = 0.;
    }

    norm1 = 1./norm1;
    norm2 = 1./norm2;
    ptf_1=0.;
    for(k_1=1;k_1<=n_1_pt;k_1++){
        foutalicept<<ptf_1+(dptf_1/2)<<","<<Yridge[0][k_1-1]*norm1<<","<<Yridge[1][k_1-1]*norm2<<","<<sum_1j<<std::endl;
        ptf_1 += dptf_1;
    }
    foutalicept.close();
    delete [] Yridge;
}


//pt1~2, 2~3, 3~4 부터 시작해보자.

void func2(double ptf_st, double ptf_end, double etaf_st, double etaf_end, double etacms_st, double etacms_end, double phif_st, double phif_end, int n, int check2){
    
    double ptf, etaf, phif, sum_alice, sum_cms, sum_j, etacms;
    double dptf, dphif, detaf, detacms, delta_Deltaeta, delta_Deltaetacms, ptf0, etaf0, phif0, etacms0;
    int i, j, k;

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
                
                sum_alice += frnk(ptf)*ptf*RidgeDis(Aridge, ptf, etaf, phif, check2)*dptf*detaf/delta_Deltaeta;
                sum_cms += frnk(ptf)*ptf*RidgeDis(Aridge, ptf, etacms, phif, check2)*dptf*detacms/delta_Deltaetacms;
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


// int n_7 = 100;

// double dptf_7 = double ((40.-10.)/n_7);
// double ptf0_7 = 10.;
// double detaf_7 = double ((1.8-1.6)/n_7);
// double etaf0_7 = 1.6;
// double dphif_7 = double ((1.28+1.28)/n_7);  //func2 phi 출력 범위 : -1.28~1.28
// double phif_7 = -1.28;
// int check2_7 = 1;




// void func7(){
//     double ptf_7, etaf_7, phif0_7, sum_7j, sum_7;
//     int i_7, j_7, k_7;
//     std::ofstream fout0002("Jet_PhiCorrelation_pT_over10.csv");
//     fout0002<<"phi,d^2N/dphi\n";
//     // double dptf;


//     for(k_7=1;k_7<=n_7+1;k_7++){
//         etaf_7 = etaf0_7;
//         // cout<<phif<<endl;
//         for(i_7=1;i_7<=n_7+1;i_7++){
//             ptf_7 = ptf0_7;
//             // cout<<etaf<<endl;
//             for(j_7=1;j_7<=n_7+1;j_7++){
//                 // sum_7 = ptf_7*RidgeDis(Aridge, ptf_7, etaf_7, phif_7, check2_7)*1/3+fj*integralNjet(ptf_7, etaf_7, phif_7, constant);
//                 // sum_7 = frnk(ptf_7)*ptf_7*fj*integralNjet(ptf_7, etaf_7, phif_7, constant)*dphif_7*detaf_7;
//                 // sum_7 *= frnk(ptf_7);
//                 // sum_7 += ptf_7*RidgeDis(Aridge, ptf_7, etaf_7, phif_7, check2_7)*dphif_7*detaf_7;
//                 // cout<<ptf_7<<setw(15)<<etaf_7<<setw(15)<<phif_7<<setw(15)<<sum_7<<endl;
//                 sum_7j += ptf_7*integralNjet(ptf_7,etaf_7,phif_7,constant)*dptf_7*detaf_7/delta_Deltaeta;
//                 ptf_7 += dptf_7;
                
//             }
            
//             etaf_7 += detaf_7;
//         }
//         sum_7 *= 2*2/3;
//         sum_7j *= 2;
//         // fprintf(fpt,"%f, %f\n", phif, sum);

//         fout0002<<phif_7<<","<<sum_7j<<std::endl;
//         sum_7 = sum_7j = 0.;
//         phif_7 += dphif_7;
//     }
//     fout0002.close();
// }


// void func8(){
//     double ptf_8, etaf_8, phif_8, detaf_8, dphif_8, etaf0_8, phif0_8, dist_8, dptf_8, sum_8, sum_8j;
//     int netaf_8, nphif_8, n_8, check2_8, i_8, j_8;
//     std::ofstream fout1012("Y^near.csv");
//     fout1012<<"pt^Jet,Y^near\n";
//     // double dptf;

//     n_8 = 100;
//     // dptf_7 = double ((1.0-0.1)/n_7);
//     detaf_8 = double ((1.2-1.0)/n_8);
//     dphif_8 = double ((M_PI+M_PI)/n_8);
//     sum_8 = 0.;
//     etaf_8 = 1.;
//     phif_8 = -M_PI;

//     check2_8 = 1;

//     for(ptf_8=0.1;ptf_8<=2.;ptf_8+=0.01){
//         etaf_8 = 1.;
//         // cout<<phif<<endl;
//         for(i_8=1;i_8<=n_8+1;i_8++){
//             phif_8 = -M_PI;
//             // cout<<etaf<<endl;
//             for(j_8=1;j_8<=n_8+1;j_8++){
//                 // sum_8 = ptf_8*RidgeDis(Aridge, ptf_8, etaf_8, phif_8, check2_8)*1/3+fj*integralNjet(ptf_8, etaf_8, phif_8, constant);
//                 // sum_8 = frnk(ptf_8)*ptf_8*fj*integralNjet(ptf_8, etaf_8, phif_8, constant)*dphif_8*detaf_8;
//                 sum_8 += ptf_8*RidgeDis(Aridge, ptf_8, etaf_8, phif_8, check2_8)*dphif_8*detaf_8;
//                 sum_8j += ptf_8*integralNjet(ptf_8,etaf_8,phif_8,constant)*dphif_8*dptf_8;
//                 // cout<<ptf_8<<setw(15)<<etaf_8<<setw(15)<<phif_8<<setw(15)<<sum_8<<endl;
//                 phif_8 += dphif_8;
//             }
//             etaf_8 += detaf_8;
//         }
//         sum_8 *= 2*2/3;
//         // fprintf(fpt,"%f, %f\n", phif, sum);

//         fout1012<<ptf_8<<","<<sum_8<<","<<sum_8j<<std::endl;
//         sum_8 = sum_8j = 0.;
//     }
//     fout1012.close();
// }



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

    // foutoh.close();
    
    // foutjet.close();
    

    thread t1(func1);   //pt distribution - alice



    thread t2(func2, 1., 2., 1.6, 1.8, 2., 4., -1.28, 1.28, 100, 1);   //1D eta correlation, pt = 1~2
    thread t3(func2, 2., 3., 1.6, 1.8, 2., 4., -1.28, 1.28, 100, 1);   //1D eta correlation, pt = 2~3
    thread t4(func2, 3., 4., 1.6, 1.8, 2., 4., -1.28, 1.28, 100, 1);   //1D eta correlation, pt = 3~4
    thread t5(func2, 1., 4., 1.6, 1.8, 2., 4., -1.28, 1.28, 100, 1);   //1D eta correlation, pt = 1~4
    // thread t7(func7);
    // thread t6(func6);
    
    // thread t8(func8);
    
    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    // t7.join();

    // t6.join();
    
    // t8.join();

    char ch[100];
    std::string command = "python3 Graphs.py";
    strcpy(ch,command.c_str());
    std::system(ch);

    // end = clock();
    end = time(NULL);
    cout<<double(end-start)<<endl;

    return 0;
}
