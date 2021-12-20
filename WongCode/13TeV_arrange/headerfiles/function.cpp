#include "function.hpp"

double Jet::integralNjet(double pt, double eta, double phi){
    double sigmaphi;
    double constant = Njet/(Tjet*(m+Tjet)*2*M_PI);
    sigmaphi = (sigmaphizero*sigmaphizero*ma*ma)/(ma*ma+pt*pt);
    return (constant/sigmaphi)*exp(((m-sqrt(m*m+pt*pt))/Tjet)-((phi*phi+eta*eta)/(2*sigmaphi)));
}

double Aridge::rapidityintit(double pt, double eta){
    double root=sqrt((pt*pt*cosh(eta)*cosh(eta))+mb*mb);
    double a = root+pt*sinh(eta);
    double b = root-pt*sinh(eta);
    return (log(a/b))/2;
}

double Aridge::lightcone(double pti, double yi){
    
    double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mp

    double squareroot=sqrt(m*m+pti*pti);
    double yiabs = std::fabs(yi);
    return (squareroot/m)*exp(yiabs-yb);
}

double Aridge::integralAridge(double pti, double yi){
    double x = lightcone(pti, yi);
    double squareroot=sqrt(m*m+pti*pti);
    if(x>=1.){
        return 0.;
    }
    else{
        return pti*pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
    }
    
}

double Ridge::lightcone(double pti, double yi){    
    double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mp

    double squareroot=sqrt(m*m+pti*pti);
    double yiabs = std::fabs(yi);
    return (squareroot/m)*exp(yiabs-yb);
}

double Ridge::RidgeDisi(double pti, double yi){
    double x = lightcone(pti, yi);
    double squareroot=sqrt(m*m+pti*pti);


        if(x>=1.){
            return 0.;
        }
        else{
            return pow(1-x,a)*exp(-sqrt(m*m+pti*pti)/T)/sqrt(md*md+pti*pti);
        }
    
}


double Ridge::RidgeDisf(double ptf, double etaf, double phif){
    double etajet = 0.;
    double ptisq = ptf*ptf-2*ptf*q*cos(phif)/cosh(etajet)+q*q/(cosh(etajet)*cosh(etajet));
    double pti;
    if(ptisq<0.0000000001){
        pti = 0.;
    }
    pti = sqrt(ptisq);
    double E = sqrt(ptf*ptf*cosh(etaf)*cosh(etaf)+m*m);
    double Ei = sqrt(pti*pti+ptf*ptf*sinh(etaf)*sinh(etaf)+m*m);

    double yi = log((Ei+ptf*sinh(etaf))/(Ei-ptf*sinh(etaf)))/2;
    double yf = log((E+ptf*sinh(etaf))/(E-ptf*sinh(etaf)))/2;

    double x = lightcone(pti, yi);

    if (x>=1.){
        return 0.;
    }
    
    else{              
        return (Aridge*RidgeDisi(pti, yi))*sqrt(1.-((mb*mb)/((mb*mb+ptf*ptf)*cosh(yf)*cosh(yf))))*(E/Ei);
    }
    
}