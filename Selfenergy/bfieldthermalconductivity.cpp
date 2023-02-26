#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include "parameter.h"

#define gbscatteringrate 0.01

using namespace std;

double dispersion(double,double,double);
double A(double,double,double);
double B(double,double,double);
double gamma(double,double);
double heatconductivitydensity(double,double,double,double,double);

int main(){
    double T = 0.1;
    int limit = (ksize)*(dksize+1);
    int kspacesize = limit-dksize;
    double* help = new double[kspacesize*kspacesize];
    while(T<=Tbfield){
    double b = 0.0;
    ostringstream fout;
    fout << "Bfieldthermalconductivity/bfieldthermalconductivityT"<< T*100 << "cJ.dat";
    ofstream write(fout.str().c_str());
    while(b<1.0){
        if(b>=compactBstart){
            ostringstream fin;
            fin << "Compactinterpolateddecayrate/interpolateddecayrateB" << b*1000 << "mBsT" << T*100 << "cJ.dat";
            ifstream read(fin.str().c_str(),ios_base::binary);
            for(int kx=0; kx<kspacesize; kx++){
                for(int ky=0; ky<kspacesize; ky++){
                    read.seekg((kx*kspacesize+ky)*sizeof(double));
                    read.read((char*)& help[kx*kspacesize+ky], sizeof(double));
                }
            }
            read.close();
            double sumnodecay = 0.0;
            double sumdecay = 0.0;
            for(int kx=0; kx<kspacesize; kx++){
                for(int ky=0; ky<kspacesize; ky++){
                    if((kx!=kspacesize-1) || (ky!=kspacesize-1)){
                        if(help[kx*kspacesize+ky]>10*gbscatteringrate){
                            sumdecay = sumdecay + heatconductivitydensity(kx*M_PI/(kspacesize-1),ky*M_PI/(kspacesize-1),b,1.0/T,J*gbscatteringrate+help[kx*kspacesize+ky])/(kspacesize*kspacesize);
                            sumnodecay = sumnodecay + heatconductivitydensity(kx*M_PI/(kspacesize-1),ky*M_PI/(kspacesize-1),b,1.0/T,J*gbscatteringrate)/(kspacesize*kspacesize);
                        }
                        else{
                            sumdecay = sumdecay + heatconductivitydensity(kx*M_PI/(kspacesize-1),ky*M_PI/(kspacesize-1),b,1.0/T,J*gbscatteringrate)/(kspacesize*kspacesize);
                            sumnodecay = sumnodecay + heatconductivitydensity(kx*M_PI/(kspacesize-1),ky*M_PI/(kspacesize-1),b,1.0/T,J*gbscatteringrate)/(kspacesize*kspacesize);
                        }
                        
                    }
                }
            }
            write << b << '\t' << sumdecay << '\t' << sumnodecay << '\n';
        }
        else{
            double sum = 0;
            for(int kx=0; kx<kspacesize; kx++){
                for(int ky=0; ky<kspacesize; ky++){
                    if((kx==kspacesize-1) || (ky==kspacesize-1)){}
                    else if((kx==0) || (ky==0)){}
                    else{
                        sum = sum + heatconductivitydensity(kx*M_PI/(kspacesize-1),ky*M_PI/(kspacesize-1),b,1.0/T,J*gbscatteringrate)/(kspacesize*kspacesize);
                    }
                }
            }
            write << b << '\t' << sum << '\n';
        }
        b = b + compactBstep;
    }
    T = T + 0.1;
    cerr << (T-Tstart)/Tend << "\r";
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

double heatconductivitydensity(double kx, double ky, double b, double beta, double decayrate){
    return((beta*beta*J*J*J*(sin(kx)*sin(kx)+sin(ky)*sin(ky))*(b*b+(2*b*b-1)*gamma(kx,ky))*(b*b+(2*b*b-1)*gamma(kx,ky)))/((decayrate)*(cosh(dispersion(kx,ky,b)*beta)-1.0)));
}
        
