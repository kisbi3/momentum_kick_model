#ifndef function_h
#define function_h
#include <cmath>




class Jet
{
    public:
        // double m;
        // double Njet;
        // double Tjet;
        // double ma;
        // double sigmaphizero;
        double integralNjet(double pt, double eta, double phi);
};


class Aridge
{
    public:
        // double mp;
        // double sqrSnn;
        // double m;
        // double mb;
        // double a;
        // double T;
        // double md;
        double rapidityintit(double pt, double eta);
        double lightcone(double pti, double yi);
        double integralAridge(double pti, double yi);
};

class Ridge
{
    public:
        // double mp;
        // double sqrSnn;
        // double m;
        // double mb;
        // double a;
        // double T;
        // double md;
        // double q;
        double Aridge;
        double lightcone(double pti, double yi);
        double RidgeDisi(double pti, double yi);
        double RidgeDisf(double ptf, double etaf, double phif);
};




#endif