#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <complex>
#include "parameter.h"

#define kxstart 0
#define kystart 0

using namespace std;

int main(){
    int count, kx, ky;
    double b = bstart;
    complex<double>* help = new complex<double>[4*(ksize)-3];
    double klength;
    while(b<1.0){
        ostringstream fin;
        fin << "selfenergy/selfenergyB" << b*1000 << "mBs.dat";
        ifstream read(fin.str().c_str(),ios_base::binary);
        count = 0;
        kx = kxstart;
        ky = kystart;
        while(count<4*(ksize)-3){
            if(count<(ksize)){
                read.seekg(((kx*(ksize)+ky))*sizeof(complex<double>));
                read.read((char*)& help[count], sizeof(complex<double>));
                if(kx<ksize-1+kxstart){
                    kx++;
                    ky++;
                }
                    count++;
            }
            else if(count>=(ksize) && count<2*(ksize)-1){
                if(kx>kxstart){
                    kx--;
                }
                read.seekg(((kx*(ksize)+ky))*sizeof(complex<double>));
                read.read((char*)& help[count], sizeof(complex<double>));
                count++;
            }
            else if(count>=2*(ksize)-1 && count <= 3*(ksize)-3){
                if(kx<ksize-1+kxstart){
                    kx++;
                    ky--;
                }
                read.seekg(((kx*(ksize)+ky))*sizeof(complex<double>));
                read.read((char*)& help[count], sizeof(complex<double>));
                count++;
            }
            else if(count>3*(ksize)-3){
                if(kx>kxstart){
                    kx--;
                }
                read.seekg(((kx*(ksize)+ky))*sizeof(complex<double>));
                read.read((char*)& help[count], sizeof(complex<double>));
                count++;
            }
        }
        read.close();
        ostringstream fout;
        fout << "selfenergyplot/decayrateplotB" << b*1000 << "mBs.dat";
        ofstream write(fout.str().c_str());
        klength = 0.0;
        for(count=0; count<4*(ksize)-3;count++){
            if(count<(ksize)){
                write << klength << '\t' << -imag(help[count]) << '\n';
                klength = klength + sqrt(2.0)*M_PI/(ksize-1);
            }
            else if(count>=(ksize) && count<2*(ksize)-1){
                write << klength << '\t' << -imag(help[count]) << '\n';
                klength = klength + M_PI/(ksize-1);
            }
            else if(count>=2*(ksize)-1 && count <= 3*(ksize)-3){
                write << klength << '\t' << -imag(help[count]) << '\n';
                klength = klength + sqrt(2.0)*M_PI/(ksize-1);
            }
            else if(count>3*(ksize)-3){
                write << klength << '\t' << -imag(help[count]) << '\n';
                klength = klength + M_PI/(ksize-1);
            }
        }
        write.close();
        ostringstream fout2;
        fout2 << "selfenergyplot/energycorrectionplotB" << b*100 << "cBs.dat";
        ofstream write2(fout2.str().c_str());
        klength = 0.0;
        for(count=0; count<4*(ksize)-3;count++){
            if(count<(ksize)){
                write2 << klength << '\t' << real(help[count]) << '\n';
                klength = klength + sqrt(2.0)*M_PI/(ksize-1);
            }
            else if(count>=(ksize) && count<2*(ksize)-1){
                write2 << klength << '\t' << real(help[count]) << '\n';
                klength = klength + M_PI/(ksize-1);
            }
            else if(count>=2*(ksize)-1 && count <= 3*(ksize)-2){
                write2 << klength << '\t' << real(help[count]) << '\n';
                klength = klength + sqrt(2.0)*M_PI/(ksize-1);
            }
            else if(count>=3*(ksize)-2){
                write2 << klength << '\t' << real(help[count]) << '\n';
                klength = klength + M_PI/(ksize-1);
            }
        }
        write2.close();
        if(b<0.7){
            b = b + 0.2;
        }
        else{
            b = b + 0.025;
        }
    }
    return 0;
}
