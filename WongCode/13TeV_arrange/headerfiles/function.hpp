#ifndef function_h
#define function_h
#include <cmath>


//static이 없을 경우 오류를 일으킴
//상수 모음

//Ridge Parameters
static double a = .5;    //fall off parameter
static double T = .62;    //Temperature, GeV
static double q = .87;    //GeV
static double m = 0.13957018;  //m == mpi
static double mb = m; //mb==mpi, GeV
static double md = 1.;   //GeV
static double sqrSnn = 13000.;
static double mp = 0.938272046; //Proton mass, GeV


//Jet Parameters
static double Njet=10.;
static double fj=1.;
static double Tjet=0.7;
static double sigmaphizero = .05;
static double ma = 100.;




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