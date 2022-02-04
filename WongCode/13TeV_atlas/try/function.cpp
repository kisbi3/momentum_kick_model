#define _USE_MATH_DEFINES

//컴파일 명령어 : g++ -shared -fPIC -o <동적 라이브러리 이름>.so ./<컴파일 대상 파일 이름>.cpp
//g++ -shared -fPIC -o function.so ./function.cpp

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
double T = .29;    //Temperature, GeV
double q_1 = .2;    //GeV
double m = 0.13957018;  //m == mpi
double mb = m; //mb==mpi, GeV
double md = 1.;   //GeV
double sqrSnn = 13000.;
double mp = 0.938272046; //Proton mass, GeV


double xx = .67;
double yy = .71;


double etajet = 0.;

double Aridge;

double frnk(double pt){
    return  xx*exp(yy*pt);
}

double rapidityintit(double pt, double eta){
    double root=sqrt((pt*pt*cosh(eta)*cosh(eta))+mb*mb);
    double a = root+pt*sinh(eta);
    double b = root-pt*sinh(eta);
    return (log(a/b))/2;
}

double lightcone(double pti, double yi){
    double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mp
    double squareroot=sqrt(m*m+pti*pti);
    double yiabs = std::fabs(yi);
    return (squareroot/m)*exp(yiabs-yb);
}


//Aridge를 구하기 위해 적분할 함수
double integralAridge(double pti, double yi, int check){ 
    double x = lightcone(pti, yi);
    double squareroot=sqrt(m*m+pti*pti);

    if(check == 1){
        if(x>=1.){
            return 0.;
        }
        else{
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
    double E = sqrt(ptf*ptf*cosh(etaf)*cosh(etaf)+m*m);
    double Ei = sqrt(pti*pti+ptf*ptf*sinh(etaf)*sinh(etaf)+m*m);

    double yi = log((Ei+ptf*sinh(etaf))/(Ei-ptf*sinh(etaf)))/2;
    double yf = log((E+ptf*sinh(etaf))/(E-ptf*sinh(etaf)))/2;

    double x = lightcone(pti, yi);

    if (x>=1.){
        return 0.;
    }
    
    else{              
        return (Aridge*integralAridge(pti, yi, check))*sqrt(1.-((mb*mb)/((mb*mb+ptf*ptf)*cosh(yf)*cosh(yf))))*(E/Ei);           // E/Ei 임을 명심하자.
    }//
    
}



// double arr[3]= {0.};

//multithread 검색했을 때에 void가 아닌 경우가 없었음 -> return 안되는듯?
//이번 경우에는 계속 더하기만 하면 되기 때문에 mutex로 더할 필요 없음. 단, 다른 메모리에서 접근하여 가져갈 위험이 있는 경우, mutex의 lock/unlock 필수!
extern "C"{
    void freeptr(void *ptr)
    {
        free(ptr);
    }
}

//pt1~2, 2~3, 3~4 부터 시작해보자.
extern "C"{
    double* py_phiCorr(double Aridge, double ptf_st, double ptf_end, double etaf_st, double etaf_end, double etacms_st, double etacms_end, double etaatlas_st, double etaatlas_end, double phif_st, double phif_end, int n, int check2, double ptjetcut, double q){
        double ptf, etaf, phif, sum_alice, sum_cms, sum_j, etacms, etaatlas, sum_atlas, norm;
        double dptf, dphif, detaf, detacms, detaatlas, delta_Deltaeta, delta_Deltaetacms, delta_Deltaetaatlas, ptf0, etaf0, phif0, etacms0, etaatlas0;
        int i, j, k;
        norm = 0.;
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

        double *result = (double*)malloc(n*sizeof(double));
        for(k=1;k<=n+1;k++){
            etacms = etacms0;
            etaf = etaf0;
            etaatlas = etaatlas0;
            for(i=1;i<=n+1;i++){
                ptf = ptf0;
                for(j=1;j<=n+1;j+=1){
                    sum_atlas += frnk(ptf)*ptf*RidgeDis(Aridge, ptf, etaatlas, phif, check2, q)*dptf*detaatlas/delta_Deltaetaatlas;
                    ptf += dptf;
                }
                etaf += detaf;
                etacms += detacms;
                etaatlas += detaatlas;
            }
            sum_atlas *= 2*2/3;

            norm += sum_atlas*dphif;
            result[k-1] = sum_atlas;

            phif += dphif;
            sum_alice = sum_cms = sum_atlas = sum_j = 0.;
        }
        if(phif <= 1.1){
            std::cout<<"phif : "<<phif<<"\tIntegrate result\t"<<norm<<std::endl;
        }

        return result;

    }
}

extern "C"{
    double py_Aridge(){
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
                sum += intaridge;

                if(i == nyi){
                    sum += (integralAridge(pti, yi0, check2)-intaridge)/2.;
                    sum *= dyi;
                    break;
                }
                else if(i != 1 && intaridge == 0.){
                    sum += (integralAridge(pti, yi0, check2)-check)/2.;
                    sum *= dyi;
                    break;       
                }

            }

            yi = 0.;
            sum *= 2.;
            double checksum;
            totalsum += sum;
            if (k==1){
                checksum = totalsum;
            }

            else if (k==npti+1 || lightcone(pti+dpti,0.)>=1.){
                totalsum -= (sum+checksum)/2.;
                totalsum *= dpti;
                break;                
            }
            pti += dpti;
            sum = 0.;
        }
        totalsum *= 2*M_PI;
        Aridge = 1/totalsum;
        return Aridge;
    }
}



// extern "C"{
//     double phicorr()
//     {
//         //multithread 컴파일 할 경우 다음의 명령어를 입력하자.
//         // g++ -o Ridge Ridge.cpp -lpthread
//         using std::cout;
//         using std::endl;
//         using std::setw;
//         using std::thread;




//         return func2(1., 2., 1.6, 1.8, 2., 4., 2., 5., -1.28, 1.28, 300, 1, 0., q_1);
//     }
// }