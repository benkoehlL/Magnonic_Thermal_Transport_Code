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

int main(){
    int kpi = ksize - 1;
    int qpi = qsize - 1;
    double* help = new double[8*qsize*qsize];
    double b = bstart;
    complex<double> sum;
    double energydiff;
    double energysum;
    complex<double> selfenergy;
    complex<double> epsilon = {0.0,delta};
    while(b<bmax){
        ostringstream fin;
        fin << "vertex/vertexB" << b*1000 << "mBs.dat";
        ifstream read(fin.str().c_str(),ios_base::binary);
        ostringstream fout;
        fout << "selfenergy/selfenergyB" << b*1000 << "mBs.dat";
        ofstream write(fout.str().c_str(),ios_base::binary);
        for(int kx=0; kx<ksize; kx++){
            for(int ky=0; ky<ksize; ky++){
                sum = {0.0,0.0};
                for(int qx=0; qx<2*qsize; qx++){
                    for(int qy=0; qy<2*qsize; qy++){
                            read.seekg((2*(((2*qsize)*(2*qsize))*(kx*(ksize)+ky)+qx*(2*qsize)+qy))*sizeof(double));
                            read.read((char*)& help[2*(qx*(2*qsize)+qy)], sizeof(double));
                            read.seekg((2*(((2*qsize)*(2*qsize))*(kx*(ksize)+ky)+qx*(2*qsize)+qy)+1)*sizeof(double));
                            read.read((char*)& help[(2*(qx*(2*qsize)+qy)+1)], sizeof(double));
                    }
                }
                for(int qx=0; qx<2*qsize; qx++){
                        for(int qy=0; qy<2*qsize; qy++){
                            energydiff = dispersion(kx*M_PI/kpi,ky*M_PI/kpi,b) - dispersion(qx*M_PI/qpi,qy*M_PI/qpi,b) - dispersion((kx*(M_PI/kpi)-qx*(M_PI/qpi))+M_PI,(ky*(M_PI/kpi)-qy*(M_PI/qpi))+M_PI,b);
                            energysum = dispersion(kx*M_PI/kpi,ky*M_PI/kpi,b) + dispersion(qx*M_PI/qpi,qy*M_PI/qpi,b) + dispersion((kx*(M_PI/kpi)+qx*(M_PI/qpi))-M_PI,(ky*(M_PI/kpi)+qy*(M_PI/qpi))-M_PI,b);
                            sum = sum + help[2*(qx*(2*qsize)+qy)]/(energydiff + epsilon);
                            sum = sum - help[2*(qx*(2*qsize)+qy)+1]/(energysum - epsilon);                            
                        }
                }
                selfenergy = sum/((complex<double>)((2*qsize)*(2*qsize)));
                write.write((char*)& selfenergy, sizeof(complex<double>));
            }
        }
        write.close();
        read.close();
        if(b<0.7){
            b = b + 0.2;
        }
        else{
            b = b + 0.025;
        }
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
