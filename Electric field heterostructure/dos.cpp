#include <iostream>
#include <cmath>
#include <complex>
#include <sstream>
#include <fstream>
#include "parameter.h"

using namespace std;

double bose(double,double);

int main(){
    double T = 0.2;
    double help;
    double* energy = new double[2*systemsize];
    complex<double> sum;
    complex<double> partition;
    complex<double> cutoffepsilon;
    cutoffepsilon = {0.0,cutoffdeltados};
    for(int D=1;D<=100;D++){
        //for(int field=0;field<=xsize;field= field+2){
        int field = 6;
            ostringstream fin;
            fin << "Bogoliubov/eigenvalueD" << D << "cJxsize" << xsize << "ysize" << ysize << "fieldsize" << field << ".dat";
            ifstream read(fin.str().c_str(),ios_base::binary);
            for(int x=0;x<2*systemsize;x++){
                read.seekg((x)*sizeof(double));
                read.read((char*)& help, sizeof(double));                     
                energy[x] = abs(help);
            }
            ostringstream of1;
            of1 << "DOS/xsize" << xsize << "ysize" << ysize << "/fieldsize" << field << "/D" << D << "cJ.dat";
            ofstream writen(of1.str().c_str());//add when using binary format: ,ios_base::binary
            
            ostringstream of2;
            of2 << "DOS/xsize" << xsize << "ysize" << ysize << "/T" << T*100 << "cJ/fieldsize" << field << "/D" << D << "cJmod.dat";
            ofstream writem(of2.str().c_str());//add when using binary format: ,ios_base::binary
            for(double omega=0.0; omega < 8.0*J*S; omega=omega+8.0*J*S/250){
                sum = 0.0;
                for(int x=0;x<2*systemsize;x++){
                    sum = sum-1.0/((omega-energy[x])+cutoffepsilon);
                }
                writen << omega << '\t' << imag(sum)/(2*systemsize) << '\n';
                writem << omega << '\t' << imag(sum)/(2*systemsize)*bose(omega,T) << '\n';
                if(omega>0.06 && omega<0.08){
                    cout << D*0.01 << '\t' << imag(sum)/(2*systemsize)*bose(omega,T) << '\n';
                }
            }
            writen.close();
            writem.close();
        //}
    }
    return 0;
}

double bose(double energy, double T){
    if(energy>=epsilon){
        return (1.0/(exp(energy/T)-1));
    }
    else{
        return 0.0;
    }
}
