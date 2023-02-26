#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "parameter.h"

using namespace std;

double u(double, double, double);
double v(double, double, double);
double gamma(double,double);
double A(double,double,double);
double B(double,double,double);
double Phi1(double,double,double,double,double);
double Phi2(double,double,double,double,double);


int main(){
    double b = bstart;
    double vertex1;
    double vertexsquare1;
    double vertex2;
    double vertexsquare2;
    int kpi = ksize - 1;
    int qpi = qsize - 1;
    while(b<bmax){
        ostringstream fn;
        fn << "vertex/vertexB" << b*1000 << "mBs.dat";
        ofstream write(fn.str().c_str(),ios_base::binary);
        for(int kx=0; kx<ksize; kx++){
            for(int ky=0; ky<ksize; ky++){
                for(int qx=0; qx<2*qsize; qx++){
                    for(int qy=0; qy<2*qsize; qy++){
                        if(((abs(kx*qpi-qx*kpi)!=2*kpi*qpi) && (abs(ky*qpi-qy*kpi)!=2*qpi*kpi) && (abs(kx*qpi-qx*kpi)!=0) && (abs(ky*qpi-qy*kpi)!=0))  && ((qx != qpi) || (qy != qpi)) && ((kx != kpi) || (ky != kpi))){
                            vertex1 = Phi1(kx*M_PI/kpi,ky*M_PI/kpi,qx*M_PI/qpi,qy*M_PI/qpi,b);
                        }
                        else{
                            vertex1=0.0;
                        }
                        if(((kx*qpi+qx*kpi!=2*qpi*kpi) && (ky*qpi+qy*kpi!=2*qpi*kpi) && (kx*qpi+qx*kpi!=0) && (ky*qpi+qy*kpi!=0) && (kx*qpi+qx*kpi!=4*qpi*kpi) && (ky*qpi+qy*kpi!=4*qpi*kpi)) && ((qx != qpi) || (qy != qpi)) && ((kx != kpi) || (ky != kpi))){
                            vertex2 = Phi2(kx*M_PI/kpi,ky*M_PI/kpi,qx*M_PI/qpi,qy*M_PI/qpi,b);
                        }
                        else{
                            vertex2=0.0;
                        }
                        vertexsquare1 = 0.5*vertex1*vertex1;
                        vertexsquare2 = 0.5*vertex2*vertex2;
                        write.write((char*)& vertexsquare1, sizeof(double));
                        write.write((char*)& vertexsquare2, sizeof(double));
                    }
                }
            }
        }
        write.close();
        if(b<0.7){
            b = b + 0.2;
        }
        else{
            b = b + 0.025;
        }
    }
    return 0;
}

double u(double kx, double ky, double b){
    return (cosh(0.5*atanh(B(kx,ky,b)/(A(kx,ky,b)))));
}
    
double v(double kx, double ky, double b){
    return (sinh(0.5*atanh(B(kx,ky,b)/(A(kx,ky,b)))));
}

double A(double kx, double ky, double b){
    return (4*J*S*(1.0 + gamma(kx,ky)*b*b));
}

double B(double kx, double ky, double b){
    return (4*J*S*(gamma(kx,ky)*(1.0-b*b)));
}

double gamma(double kx, double ky){
    return (0.5*(cos(kx)+cos(ky)));
}

double dispersion(double kx, double ky, double b){
    double a = A(kx,ky,b);
    double c = B(kx,ky,b);
    return(sqrt(a*a-c*c));
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
