#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#define step 0.02

using namespace std;

double u(double, double, double);
double v(double, double, double);
double gamma(double,double);
double omega(double,double,double);
double A(double,double,double);
double B(double,double,double);

int main(){
    double b =0.5;
    fstream utest("utest.dat",ios::out);
     if (utest.good())
    {
        cout << "Alles gut gegangen..." << '\n';
    }
    else {
        cout << "Oh je...";
        utest.clear();
    }
    for(double x=0.0; x<M_PI; x=x+step*M_PI){
        utest << x << '\t' << u(x,x,b) << '\n';
    }
    utest.close();
    fstream vtest("vtest.dat",ios::out);
    for(double x=0.0; x<M_PI; x=x+step*M_PI){
        vtest << x << '\t' << v(x,x,b) << '\n';
    }
    vtest.close();
    return 0;
}


double u(double kx, double ky, double b){
    return (cosh(atanh(B(kx,ky,b)/(A(kx,ky,b)))*0.5));
}

double v(double kx, double ky, double b){
    return (sinh(atanh(B(kx,ky,b)/(A(kx,ky,b)))*0.5));
}

double A(double kx, double ky, double b){
    return (1.0 + gamma(kx,ky)*b*b);
}

double B(double kx, double ky, double b){
    return (gamma(kx,ky)*(1-b*b));
}

double gamma(double kx, double ky){
    return (0.5*(cos(kx)+cos(ky)));
}

double omega(double kx,double ky,double b){
    return (sqrt(A(kx,ky,b)*A(kx,ky,b)-B(kx,ky,b)*B(kx,ky,b)));
}
