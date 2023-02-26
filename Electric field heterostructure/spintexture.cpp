#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include "parameter.h"

using namespace std;

int main(){
    
    double d = 0.01;
    double* state = new double[systemsize];
    while(d<1.0){
    ostringstream fin;
    fin << "Texture/textureD" << d*100 << "cJxsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
    ifstream read(fin.str().c_str(),ios_base::binary);
    for(int x=0;x<xsize;x++){
        for(int y=0;y<ysize;y++){
            read.seekg((x*ysize+y)*sizeof(double));
            read.read((char*)& state[x*ysize+y], sizeof(double));
        }
    }
    read.close();
    ostringstream fout;
    fout << "PlotSpintexture/SpinTexture"<< d*100 << "cJxsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
    ofstream write(fout.str().c_str());
    for(int x=0;x<xsize;x++){
        for(int y=0;y<ysize;y++){
            write << state[x*ysize+y]+M_PI*((y+x)%2) << '\t';
        }
        write << '\n';
    }
    write.close();
    d += 0.01;
    }
    return 0;
}
