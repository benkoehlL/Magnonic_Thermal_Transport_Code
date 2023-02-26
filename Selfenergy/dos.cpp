#include <iostream>
#include <cmath>
#include <complex>
#include <sstream>
#include <fstream>
#include "parameter.h"

using namespace std;

double dispersion(double,double,double);
double A(double,double,double);
double B(double,double,double);
double gamma(double,double);

int main(){
    complex<double> sum;
    double b = 0.9;
    complex<double> partition;
    complex<double> epsilon = {0.0,delta};
    ofstream write("dos.dat");
    for(double omega=0.0; omega < 8.0*J*S; omega=omega+8.0*J*S/100){
        sum=0.0;
        for(int kx=0;kx<qsize; kx++){
            for(int ky=0;ky<qsize; ky++){
                sum=sum-1.0/((omega-4*J*S*dispersion(kx*M_PI/qsize,ky*M_PI/qsize,b))+epsilon);
            }
        }
        write << omega << '\t' << imag(sum)/((qsize)*(qsize)) << '\n';
    }
    return 0;
}
            

double gamma(double kx,double ky){
    return(0.5*(cos(kx)+cos(ky)));
}

double A(double kx, double ky, double b){
    return (1.0 + gamma(kx,ky)*b*b);
}

double B(double kx, double ky, double b){
    return (gamma(kx,ky)*(1.0-b*b));
}

double dispersion(double kx, double ky, double b){
    return(sqrt(A(kx,ky,b)*A(kx,ky,b)-B(kx,ky,b)*B(kx,ky,b)));
}
