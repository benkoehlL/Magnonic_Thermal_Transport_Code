#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "parameter.h"

using namespace std;

int main(){
    double b = compactBstart;
    int limit = (ksize)*(dksize+1);
    int kspacesize = limit-dksize;
    double* help = new double[kspacesize*kspacesize];
    while(b<1.0){
        ostringstream fin;
        fin << "Compactinterpolateddecayrate/interpolateddecayrateB" << b*1000 << "mBs.dat";
        ifstream read(fin.str().c_str(),ios_base::binary);
        for(int kx=0;kx<kspacesize;kx++){
            for(int ky=0;ky<kspacesize;ky++){
                read.seekg((kx*(kspacesize)+ky)*sizeof(double));
                read.read((char*)& help[kx*kspacesize+ky], sizeof(double));
            }
        }
        read.close();
        ostringstream fout;
        fout << "Densityplot/densityplotB" << b*1000 << "mBs.dat";
        ofstream write(fout.str().c_str());
        for(int kx=0;kx<kspacesize;kx++){
            for(int ky=0;ky<kspacesize;ky++){
                write << kx*M_PI/(kspacesize-1) << '\t' << ky*M_PI/(kspacesize-1) << '\t' << help[kx*kspacesize+ky] << '\n';
                //write << kx << '\t' << ky << '\t' << help[kx*kspacesize+ky] << '\n';
            }
        }
        write.close();
        /*if(b<0.7){
            b = b + 0.2;
        }
        else{
            b = b + 0.025;
        }*/
        b = b + 0.01;
    }
    return 0;
}
