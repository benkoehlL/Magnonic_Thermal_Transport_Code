#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "parameter.h"

using namespace std;

int main(){
    double *help = new double[2*systemsize];
    for(int E=0;E<=100;E++){
        ostringstream fin;
        fin << "Texture/xsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << 
        "/textureD" << E << "cJxsize_simplex.dat";
        ifstream read(fin.str().c_str(),ios_base::binary);
        for(int x=0;x<xsize;x++){
            for(int y=0;y<ysize;y++){
                read.seekg((x*ysize+y)*sizeof(double));
                read.read((char*)& help[x*ysize+y], sizeof(double));
                read.seekg((x*ysize+y+systemsize)*sizeof(double));
                read.read((char*)& help[x*ysize+y+systemsize], sizeof(double));
            }
        }
        read.close();
        ostringstream foutphi,fouttheta;
        foutphi << "Angleplot/xsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << 
        "/phiplotD" << E << "cJ.dat";
        ofstream writephi(foutphi.str().c_str());
        fouttheta << "Angleplot//xsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << 
        "/thetaplotD" << E << "cJ.dat";
        ofstream writetheta(fouttheta.str().c_str());
        for(int y=0;y<ysize;y++){
            for(int x=0;x<xsize;x++){
                while(help[x*ysize+y]<0){
                    help[x*ysize+y] = help[x*ysize+y] + 2*M_PI;
                }
                while(help[x*ysize+y+systemsize]<0){
                    help[x*ysize+y+systemsize] = help[x*ysize+y+systemsize] + 2*M_PI;
                }
                while(help[x*ysize+y]>2*M_PI){
                    help[x*ysize+y] = help[x*ysize+y] - 2*M_PI;
                }
                while(help[x*ysize+y+systemsize]>2*M_PI){
                    help[x*ysize+y+systemsize] = help[x*ysize+y+systemsize] - 2*M_PI;
                }
                /*if(help[x*ysize+y+systemsize]>M_PI){
                    help[x*ysize+y+systemsize] = 2*M_PI-help[x*ysize+y+systemsize];
                    if(help[x*ysize+y]+M_PI>2*M_PI){
                        help[x*ysize+y] = -M_PI+help[x*ysize+y];
                    }
                    else{
                        help[x*ysize+y] = +M_PI+help[x*ysize+y];
                    }
                }*/
                writephi << help[x*ysize+y] << '\t';
                writetheta << help[x*ysize+y+systemsize] << '\t';
            }
            writephi << '\n';
            writetheta << '\n';
        }
        writephi.close();
        writetheta.close();
    }
    return 0;
}
