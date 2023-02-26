#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#define size 200
#define S 0.5
#define J 1

using namespace std;

double dispersion(double,double,double);
double A(double,double,double);
double B(double,double,double);
double gamma(double,double);

int main(){
    double b = 0.0;
    double help;
    while(b<1.0){
        ostringstream fout;
        fout << "Dispersion/dispB" << b*100 << "cBs.dat";
        ofstream write(fout.str().c_str(),ios_base::binary);
        for(int kx=0; kx<size+1; kx++){
            for(int ky=0; ky<size+1; ky++){
                help = dispersion(kx*M_PI/size,ky*M_PI/size,b);
                write.write((char*)& help, sizeof(double));
            }
        }
        write.close();
        b= b + 0.1;
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
/*B-field dispersion:
double dispersion(double kx, double ky, double b){
    return(4*J*S*sqrt(A(kx,ky,b)*A(kx,ky,b)-B(kx,ky,b)*B(kx,ky,b)));
}*/
//E-field disersion
double dispersion(double kx, double ky, double b){
    return(4*J*S*sqrt((1.0-gamma(kx,ky))*(b*b+J*J+J*sqrt(b*b+J*J)*gamma(kx,ky))));
}


