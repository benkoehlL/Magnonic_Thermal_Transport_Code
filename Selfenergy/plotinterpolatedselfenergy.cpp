#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "parameter.h"

#define kxstart 0
#define kystart 0

using namespace std;

int main(){
    int count, kx, ky;
    double b = bstart;
    int limit = (ksize)*(dksize+1);
    int kspacesize = limit-dksize;
    double* help = new double[4*(kspacesize)-3];
    double klength;
    while(b<1.0){
        ostringstream fin;
        fin << "interpolatedlifetime/interpolatedlifetimeB" << b*100 << "cBs.dat";
        ifstream read(fin.str().c_str(),ios_base::binary);
        count = 0;
        kx = kxstart;
        ky = kystart;
        while(count<4*(kspacesize)-3){
            if(count<kspacesize){
                read.seekg((kx*(kspacesize)+ky)*sizeof(double));
                read.read((char*)& help[count], sizeof(double));
            if(kx<kspacesize-1+kxstart){
                kx++;
                ky++;}
                count++;
            }
            else if(count>=kspacesize && count<2*(kspacesize)-1){
                if(kx>kxstart){
                    kx--;
                }
                read.seekg((kx*(kspacesize)+ky)*sizeof(double));
                read.read((char*)& help[count], sizeof(double));
                count++;
            }
            else if(count>=2*kspacesize-1 && count <= 3*kspacesize-3){
                if(kx<kspacesize-1+kxstart){
                    kx++;
                    ky--;
                }
                read.seekg((kx*kspacesize+ky)*sizeof(double));
                read.read((char*)& help[count], sizeof(double));
                count++;
            }
            else if(count>3*kspacesize-3){
                if(kx>kxstart){
                    kx--;
                }
                read.seekg((kx*(kspacesize)+ky)*sizeof(double));
                read.read((char*)& help[count], sizeof(double));
                count++;
            }
        }
        read.close();
        ostringstream fout;
        fout << "interpolatedlifetimeplot/interpolatedlifetimeplotB" << b*100 << "cBs.dat";
        ofstream write(fout.str().c_str());
        klength = 0.0;
        for(count=0; count<4*kspacesize-3;count++){
            if(count<kspacesize+dksize){
                write << klength << '\t' << help[count] << '\n';
                klength = klength + sqrt(2.0)*M_PI/(kspacesize-1);
            }
            else if(count>=kspacesize+dksize && count<2*kspacesize+dksize-1){
                write << klength << '\t' << help[count] << '\n';
                klength = klength + M_PI/(kspacesize-1);
            }
            else if(count>=2*kspacesize+dksize-1 && count < 3*kspacesize-2+dksize){
                write << klength << '\t' << help[count] << '\n';
                klength = klength + sqrt(2.0)*M_PI/(kspacesize-1);
            }
            else if(count>=3*kspacesize-2+dksize){
                write << klength << '\t' << help[count] << '\n';
                klength = klength + M_PI/(kspacesize-1);
            }
        }
        write.close();
        b = b + 0.1;
    }
    return 0;
}
