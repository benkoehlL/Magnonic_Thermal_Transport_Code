#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
#include "parameter.h"

using namespace std;

double C(double);
double F1(double);
double F2(double);
double gamma(double,double);
double A(double,double,double);
double B(double,double,double);
double energy(double,double,double);
double bosefunction(double, double);

int main(){
    //omega behaviour without temperature dependence
    int kxpi = xsize;
    int kypi = ysize;
    complex<double> omegasumxx,omegasumxy,omegasumyx,omegasumyy;
    complex<double> cutoffepsilon;
    cutoffepsilon = {0.0,cutoffdelta};
    double prefactorx,prefactory;
    double bose;
    double omega = 0.0;
    ostringstream of1;
    of1 << "Drudeweightoffield/drudeanalyticxsize" << xsize << "ysize" << ysize << ".dat";
    ofstream writedrudeofD(of1.str().c_str());//add when using binary format: ,ios_base::binary
    for(double d=0.00;d<1.01;d+=0.01){
        ostringstream of;
        of << "Drudeweight/drudeweightanaD" << d*100 << "cJxsize" << xsize << "ysize" << ysize << ".dat";
        ofstream writedrude(of.str().c_str());//add when using binary format: ,ios_base::binary
        for(double T=0.001;T<3.0;T+=0.001){
        //for(double omega=-0.001;omega<4.1;omega+=0.001){
            omegasumxx = {0.0,0.0};
            omegasumxy = {0.0,0.0};
            omegasumyx = {0.0,0.0};
            omegasumyy = {0.0,0.0};
            for(int kx=0; kx<xsize; kx++){
                for(int ky=0; ky<ysize; ky++){
                    prefactorx = 2.0*S*S*J*J*(F1(d)+gamma(2*kx*M_PI/(kxpi),2*ky*M_PI/(kypi))*F2(d))*sin(2*kx*M_PI/(kxpi));
                    prefactory = 2.0*S*S*J*J*(F1(d)+gamma(2*kx*M_PI/(kxpi),2*ky*M_PI/(kypi))*F2(d))*sin(2*ky*M_PI/(kypi));
                    if(energy(2*kx*M_PI/(kxpi),2*ky*M_PI/(kypi),d)>=epsilon){
                        bose = bosefunction(energy(2*kx*M_PI/(kxpi),2*ky*M_PI/(kypi),d),T);
                    }
                    else{
                        bose=0.0;
                    }                
                    omegasumxx += bose*(bose+1.0)*prefactorx*prefactorx/(omega  + cutoffepsilon);
                    omegasumxy += bose*(bose+1.0)*prefactorx*prefactory/(omega  + cutoffepsilon);
                    omegasumyx += bose*(bose+1.0)*prefactory*prefactorx/(omega  + cutoffepsilon);
                    omegasumyy += bose*(bose+1.0)*prefactory*prefactory/(omega  + cutoffepsilon);
                    
                    //omegasum += prefactorxx*prefactorxx/(omega + cutoffepsilon);
            
                }
            }
            writedrude << T << '\t' << -imag(omegasumxx/((complex<double>)(systemsize))) << '\t' << -imag(omegasumxy/((complex<double>)(systemsize))) << '\t' << -imag(omegasumyx/((complex<double>)(systemsize))) << '\t' << -imag(omegasumyy/((complex<double>)(systemsize))) << '\n';

            if(T>=0.2999 && T<=0.3001){
                omegasumxx = {0.0,0.0};
                omegasumxy = {0.0,0.0};
                omegasumyx = {0.0,0.0};
                omegasumyy = {0.0,0.0};
                for(int kx=0; kx<xsize; kx++){
                    for(int ky=0; ky<ysize; ky++){
                        prefactorx = 2.0*S*S*J*J*(F1(d)+gamma(2*kx*M_PI/(kxpi),2*ky*M_PI/(kypi))*F2(d))*sin(2*kx*M_PI/(kxpi));
                        prefactory = 2.0*S*S*J*J*(F1(d)+gamma(2*kx*M_PI/(kxpi),2*ky*M_PI/(kypi))*F2(d))*sin(2*ky*M_PI/(kypi));
                        if(energy(2*kx*M_PI/(kxpi),2*ky*M_PI/(kypi),d)>=epsilon){
                            bose = bosefunction(energy(2*kx*M_PI/(kxpi),2*ky*M_PI/(kypi),d),T);
                        }
                        else{
                            bose=0.0;
                        }                
                        omegasumxx += bose*(bose+1.0)*prefactorx*prefactorx/(omega  + cutoffepsilon);
                        omegasumxy += bose*(bose+1.0)*prefactorx*prefactory/(omega  + cutoffepsilon);
                        omegasumyx += bose*(bose+1.0)*prefactory*prefactorx/(omega  + cutoffepsilon);
                        omegasumyy += bose*(bose+1.0)*prefactory*prefactory/(omega  + cutoffepsilon);
                    }
                }
                writedrudeofD << d << '\t' << -imag(omegasumxx/((complex<double>)(systemsize))) << '\t' << -imag(omegasumxy/((complex<double>)(systemsize))) << '\t' << -imag(omegasumyx/((complex<double>)(systemsize))) << '\t' << -imag(omegasumyy/((complex<double>)(systemsize))) << '\n';
            }
        }
        writedrude.close();
        cerr << d << "\r";
    }
    writedrudeofD.close();
    return 0;
}

double C(double d){
    return(sqrt(1+d*d));
}

double F1(double d){
    double c = C(d);
    return(c*(c-1));
}

double F2(double d){
    return(-2*C(d));
}

double gamma(double kx,double ky){
    return(0.5*(cos(kx)+cos(ky)));
}

double A(double kx,double ky,double d){
    double c = C(d);
    return(c+0.5*gamma(kx,ky)*(c-1.0));
}

double B(double kx,double ky, double d){
    return(-0.5*gamma(kx,ky)*(C(d)+1.0));
}

double energy(double kx,double ky,double d){
    double a = A(kx,ky,d);
    double b = B(kx,ky,d);
    return(4*S*J*sqrt(a*a-b*b));
}

double bosefunction(double energy, double T){
    return(1.0/(exp(energy/T)-1));
}
