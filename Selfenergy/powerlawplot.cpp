#include <iostream>
#include <sstream>
#include <fstream>
#include <complex>
#include "parameter.h"

using namespace std;

int main(){
    double T = 0.01;
    complex<double> help;
    ostringstream fout;
    fout << "SelfenergyB90cBsvariousTemperatures/powerlaw.dat";
    ofstream write(fout.str().c_str());
    while(T<0.8){
        ostringstream fin;
        fin << "SelfenergyB90cBsvariousTemperatures/tempselfenergyT" << T*100 << "cJ.dat";
        ifstream read(fin.str().c_str(),ios_base::binary);
        read.seekg((ksize*ksize-3)*sizeof(complex<double>));
        read.read((char*)& help, sizeof(complex<double>));
        write << T << '\t' << -imag(help) << '\n';
        T = T + 0.01;
        read.close();
    }
    write.close();
    return 0;
}
