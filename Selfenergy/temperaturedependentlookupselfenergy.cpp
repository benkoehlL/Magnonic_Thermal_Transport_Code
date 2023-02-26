#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
#include "parameter.h"

using namespace std;

double dispersion(double,double,double);
double A(double,double,double);
double B(double,double,double);
double gamma(double,double);
double u(double, double, double);
double v(double, double, double);
double Phi1(double,double,double,double,double);
double Phi2(double,double,double,double,double);
double bosefunction(double,double,double, double);

int main(){
    int kpi = ksize - 1;
    int qpi = qsize - 1;
    int temperaturecutradius = 1;
    double b = compactBstart;
    double T[3] = {0.01,0.1,0.3};
    complex<double> sum;
    double energy;
    double vertex1reversed;
    complex<double> nonvertex;
    complex<double> selfenergy;
    complex<double> epsilon = {0.0,delta};
    while(b<1.0){
        for(int i=0;i<3;i++){
            ostringstream fout;
            fout << "Compactselfenergy/temperaturedependendselfenergyB" << b*1000 << "mBsT" << T[i]*100 << "cJ.dat";
            ofstream write(fout.str().c_str(),ios_base::binary);
            for(int kx=0; kx<ksize; kx++){
                for(int ky=0; ky<ksize; ky++){
                    sum = {0.0,0.0};
                    for(int qx=0; qx<2*qsize; qx++){
                        for(int qy=0; qy<2*qsize; qy++){
                            if(((abs(kx*qpi-qx*kpi)!=2*kpi*qpi) && (abs(ky*qpi-qy*kpi)!=2*qpi*kpi) && (abs(kx*qpi-qx*kpi)!=0) && (abs(ky*qpi-qy*kpi)!=0))  && ((qx != qpi) || (qy != qpi)) && ((kx != kpi) || (ky != kpi))){
                                energy = dispersion(kx*M_PI/kpi,ky*M_PI/kpi,b) + dispersion(qx*M_PI/qpi,qy*M_PI/qpi,b) - dispersion((kx*(M_PI/kpi)-qx*(M_PI/qpi))+M_PI,(ky*(M_PI/kpi)-qy*(M_PI/qpi))+M_PI,b);
                                nonvertex = (bosefunction(qx*M_PI/qpi,qy*M_PI/qpi,b,1.0/T[i])-bosefunction((kx*(M_PI/kpi)-qx*(M_PI/qpi))+M_PI,(ky*(M_PI/kpi)-qy*(M_PI/qpi))+M_PI,b,1.0/T[i]))/(energy + epsilon);
                                vertex1reversed = Phi1(qx*M_PI/qpi,qy*M_PI/qpi,kx*M_PI/kpi,ky*M_PI/kpi,b);
                            }
                            else{
                                nonvertex = 0.0;//
                            }
                            if((qx < qsize-temperaturecutradius*double(qsize/ksize) || qx > qsize+temperaturecutradius*double(qsize/ksize)) || (qy < qsize-temperaturecutradius*double(qsize/ksize) || qy > qsize+temperaturecutradius*double(qsize/ksize))){
                                sum = sum + (0.5*vertex1reversed*vertex1reversed)*nonvertex;
                            }
                        }
                    }
                    selfenergy = sum/((complex<double>)((2*qsize)*(2*qsize)));
                    write.write((char*)& selfenergy, sizeof(complex<double>));
                }
            }
            
            write.close();
        }
        b = b + compactBstep;
    }
    return 0;
}

double gamma(double kx,double ky){
    return(0.5*(cos(kx)+cos(ky)));
}

double A(double kx, double ky, double b){
    return (4*J*S*(1.0 + gamma(kx,ky)*b*b));
}

double B(double kx, double ky, double b){
    return (4*J*S*(gamma(kx,ky)*(1.0-b*b)));
}

double dispersion(double kx, double ky, double b){
    double a = A(kx,ky,b);
    double c = B(kx,ky,b);
    return(sqrt(a*a-c*c));
} 
double u(double kx, double ky, double b){
    return (cosh(0.5*atanh(B(kx,ky,b)/(A(kx,ky,b)))));
}
    
double v(double kx, double ky, double b){
    return (sinh(0.5*atanh(B(kx,ky,b)/(A(kx,ky,b)))));
}

double Phi1(double kx, double ky, double qx, double qy, double b){
    double uk = u(kx,ky,b);
    double vk = v(kx,ky,b);
    double uq = u(qx,qy,b);
    double vq = v(qx,qy,b);
    double uQ = u(kx-qx+M_PI,ky-qy+M_PI,b);
    double vQ = v(kx-qx+M_PI,ky-qy+M_PI,b);
    double gk = gamma(kx,ky);
    double gq = gamma(qx,qy);
    double gQ = gamma(kx-qx+M_PI,ky-qy+M_PI);
    return (8*J*S*b*sqrt(1.0-b*b)/(sqrt(2.0*S))*(gk*(uk+vk)*(uq*vQ+vq*uQ)+gq*(uq+vq)*(uk*uQ+vk*vQ)+gQ*(uQ+vQ)*(uk*uq+vk*vq)));
}

double Phi2(double kx, double ky, double qx, double qy, double b){
    double uk = u(kx,ky,b);
    double vk = v(kx,ky,b);
    double uq = u(qx,qy,b);
    double vq = v(qx,qy,b);
    double uQ = u(kx+qx-M_PI,ky+qy-M_PI,b);
    double vQ = v(kx+qx-M_PI,ky+qy-M_PI,b);
    double gk = gamma(kx,ky);
    double gq = gamma(qx,qy);
    double gQ = gamma(kx+qx-M_PI,ky+qy-M_PI);
    return (8*J*S*b*sqrt(1.0-b*b)/(sqrt(2.0*S))*(gk*(uk+vk)*(uq*vQ+vq*uQ)+gq*(uq+vq)*(uk*vQ+vk*uQ)+gQ*(uQ+vQ)*(uk*vq+vk*uq)));
}

double bosefunction(double kx, double ky, double b, double beta){
    return(1.0/(exp(beta*(dispersion(kx,ky,b)))-1.0));
}
